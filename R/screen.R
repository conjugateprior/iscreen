
#' xdiscrete
#'
#' @param x
#'
#' @return
#' @export
#' @importFrom stats kmeans
#'
xdiscrete <- function(x) {
  if (length(unique(x)) > 3) {
    temp <- kmeans(x, 3, nstart = 25)
    res <- temp$cluster - 1
  } else {
    res <- x
  }

  res
}

#' makeDiscrete
#'
#' @param C
#'
#' @return Results
#' @import doParallel
#' @importFrom foreach foreach %dopar%
#' @export
makeDiscrete <- function(C) {
  cc <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(cc,
                                type = ifelse(.Platform$OS.type == "unix",
                                              "FORK", "PSOCK"))
  registerDoParallel(cl)
  cat("Discretizing variables... \n")
  res <- foreach(xx = 1:ncol(C),
                 .packages = "Matrix",
                 .export = "xdiscrete",
                 .combine = "cbind") %dopar% {
    xdiscrete(C[, xx])
  }
  colnames(res) <- colnames(C)
  parallel::stopCluster(cl)

  res
}

#' Calculate I for All Marginals
#'
#' @param C
#' @param y

#' @return Column I values, sorted from large to small
#' @import doParallel
#' @importFrom foreach foreach %dopar%
#' @export
calcI <- function(C, y) {
  res <- foreach(xx = 1:ncol(C),
                 .packages = "Matrix",
                 .export = "f.list.I",
                 .combine = c) %dopar% {
    out <- f.list.I(1, as.matrix(C[, xx]), as.vector(y))[2]
  }
  names(res) <- colnames(C)
  res <- sort(res, decreasing = TRUE)

  res
}

#' Calculate I for All Two-Way Interactions
#'
#' @param C
#' @param y
#' @param ind
#'
#' @return
#' @import doParallel
#' @importFrom foreach foreach %dopar%
#' @export
calcI2 <- function(C, y, ind) { # ind is combn(k,2) matrix
  res <- foreach(xx = 1:ncol(ind),
                 .packages = "Matrix",
                 .export = "f.list.I",
                 .combine = c) %dopar% { # rbind
    out <- f.list.I(ind[, xx], as.matrix(C), as.vector(y))[3]
  }
  r <- rbind(ind, res)
  r <- r[, which(r[3, ] >= 1)]

  r
}


#' getX2
#'
#' @param x
#'
#' @return
getX2 <- function(x) {
  name <- paste(colnames(x)[1], colnames(x)[2], sep = ":")
  X <- Matrix::Matrix(x[, 1] * x[, 2], sparse = TRUE)
  colnames(X) <- name

  X
}


#' get.predictors.int
#'
#' @param x
#' @param y
#' @param inds
#' @param nclus
#'
#' @import doParallel
#' @importFrom foreach foreach %dopar%
#' @return
get.predictors.int <- function(x, y, inds, nclus = 200) {
  cc <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(cc,
                              type = ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK"))
  registerDoParallel(cl)
  res <- foreach(j = 1:nclus,
                 .packages = "Matrix",
                 .export = "getX2",
                 .combine = cbind) %dopar% {
    getX2(x[, inds[1:2, j]])
  }
  parallel::stopCluster(cl)

  res
}

#' f.list.I
#'
#' @param var.list
#' @param data.x
#' @param data.y
#'
#' @return
#' @export
f.list.I <- function(var.list, data.x, data.y) {
  kk <- length(var.list)
  if (kk > 1) {
    xx <- as.matrix(data.x[, as.vector(var.list)] %*% as.vector((3^(0:(kk - 1)))))
  }
  else {
    xx <- as.matrix(data.x)[, var.list]
  }
  yy <- unlist(data.y)
  dat.mat <- table(xx, yy)
  n.d <- dat.mat[, 1]
  n.u <- dat.mat[, 2]
  nn.d <- sum(n.d)
  nn.u <- sum(n.u)
  i.score <- nn.d * nn.u * sum((n.d / nn.d - n.u / nn.u)^2) / (nn.d + nn.u)

  c(var.list, i.score)
}


#' f.list.screen.I
#'
#' @param var.list
#' @param data.x
#' @param data.y
#'
#' @return
#' @export
#' @importFrom utils combn
f.list.screen.I <- function(var.list, data.x, data.y) {
  kk <- length(var.list)
  mk.ind <- rep(1, kk)
  I.score <- NA
  score.pre <- f.list.I(var.list, data.x, data.y)[length(var.list) + 1]
  if (kk > 1) {
    var.list.use <- 1:kk
    data.x.use <- data.x[, var.list]
    result.v <- t(combn(var.list.use, kk - 1, f.list.I,
      simplify = T,
      data.x = data.x.use, data.y = data.y
    ))
    # print(c(score.pre, max(result.v[,kk])))
    if (max(result.v[, kk]) > score.pre) { # if after dropping each index we find one dropped version that is higher than prescore
      mk.ind[-result.v[which.max(result.v[, kk]), 1:(kk - 1)]] <- 0 # index of 0/1 turning on if list of variables kept or not
      var.list.use <- 1:(kk - 1)
      data.x.use <- data.x[, var.list[mk.ind > 0]]
      out.recur <- f.list.screen.I(var.list.use, data.x.use, data.y)
      mk.ind[mk.ind > 0] <- mk.ind[mk.ind > 0] * out.recur[(kk - 1) + (1:(kk - 1))]
      I.score <- out.recur[length(out.recur)]
    }
    else {
      I.score <- score.pre
    }
  }
  else {
    I.score <- score.pre
  }
  # if(var.list[kk]==var.eva3[length(var.eva3)]){print(c(var.list, mk.ind, I.score, date()))}
  # kept<-var.list[which(mk.ind==1)]

  c(var.list, mk.ind, I.score)
}

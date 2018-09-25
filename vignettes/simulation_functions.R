## data generation function -----------------------------------------
sim_dgp <- function(
  nobs, kcov, outcome, #treatment, 
  var_type, seed, split = FALSE
) {
  
  set.seed(seed)
  
  ## generate covariates
  if (var_type == "bern") {
    ## Bern variables
    pr <- runif(kcov, 0.3, 0.6)
    XX <- matrix(NA, nrow = nobs, ncol = kcov)
    for (k in 1:kcov) {
      XX[,k] <- rbinom(nobs, size = 1, prob = pr[k])
    }
  } else if (var_type == "normal") {
    sigma <- diag(kcov)
    sigma[sigma == 0] <- 0.3
    XX <- MASS::mvrnorm(nobs, mu = rep(0, kcov), Sigma = sigma)
  } else if (var_type == "mix") {
    kn <- ceiling(kcov/2)
    kb <- kcov - kn
    
    ## normal 
    sigma <- diag(kn)
    sigma[sigma == 0] <- 0.3
    XX1 <- MASS::mvrnorm(nobs, mu = rep(0, kn), Sigma = sigma)
    
    ## bern
    pr <- runif(kcov, 0.3, 0.6)
    XX2 <- matrix(NA, nrow = nobs, ncol = kb)
    for (k in 1:kb) { XX2[,k] <- rbinom(nobs, size = 1, prob = pr[k]) }
    
    ## combine 
    XX <- cbind(XX1, XX2)
  }
  
  ## generate outcomes
  if (outcome == "LN") {
    beta <- rnorm(kcov)
    pscore <- 1 / (1 + exp(-XX[,1:3]%*%beta[1:3]))
    influential_var<- NULL
    influential_var_indiv <- c("X1","X2","X3")
  } else if (outcome == "NL") {
    ## non-linear
    beta <- rbinom(kcov, size = 1, prob = 0.4)
    fX <- -0.3 + sin(XX[,1] * XX[,2])  - XX[,3]^2
    pscore <- 1 / (1 + exp(-fX))
    influential_var<- c("X1*X2","X3")
    influential_var_indiv <- c("X1","X2","X3")
  } else if (outcome == "LNJ") {
    ## linear with joint
    beta <- rnorm(kcov)
    pscore <- 1 / (1 + exp(-0.2-XX[,1]*XX[,2]*beta[1]) + XX[,3]*XX[,4])
    influential_var<- c("X1*X2","X3*X4")
    influential_var_indiv <- c("X1","X2","X3","X4")
  } else if (outcome == "NLJ") {
    ## non-linear with joint
    fX <- -0.1 + sin(XX[,1] * XX[,2]) - cos(XX[,4]*XX[,3])
    pscore <- 1 / (1 + exp(-fX))
    influential_var<- c("X1*X2","X3*X4")
    influential_var_indiv <- c("X1","X2","X3","X4")
  } else if (outcome == "LNMJ") {
    ## linear with marginal + joint
    beta <- rnorm(kcov)
    pscore <- 1 / (1 + exp(-0.3-XX[,3:4]%*%beta[3:4] + XX[,1]*XX[,2]))
    influential_var<- c("X1*X2","X3","X4")
    influential_var_indiv <- c("X1","X2","X3","X4")
  } else { #"NLMJ"
    ## non-linear with marginal + joint
    beta <- rbinom(kcov, size = 1, prob = 0.4)
    fX <- -0.1 + sin(XX[,1] * XX[,2]) - cos(XX[,4]*XX[,3]) + XX[,5] * beta[5]
    pscore <- 1 / (1 + exp(-fX))
    influential_var<- c("X1*X2","X3*X4","X5")
    influential_var_indiv <- c("X1","X2","X3","X4","X5")
  }
  
  Y <- rbinom(nobs, size = 1, prob = pscore)
  
  if (split) {
    half_id <- c(rep(TRUE, floor(nobs / 2)), rep(FALSE, nobs - floor(nobs / 2)))
    colnames(XX) <- paste("X", 1:kcov, sep = "")
    out_dat <- data.frame(cbind(Y,XX))
    colnames(out_dat)[1] <- c("Y")
    
    d1 <- out_dat[half_id,]
    d2 <- out_dat[!half_id,]
    
    X1 <- XX[half_id,]
    X2 <- XX[!half_id,]
    
    pscore1 <- pscore[half_id]
    pscore2 <- pscore[!half_id]
    
    out <- list(
      list("data_full" = d1, "Xmat" = X1, "pscore" = pscore1),
      list("data_full" = d2, "Xmat" = X2, "pscore" = pscore2),
      list("influential_var" = influential_var,"influential_var_indiv" = influential_var_indiv)
    )
    
  } else {
    colnames(XX) <- paste("X", 1:kcov, sep = "")
    out_dat <- data.frame(cbind(Y, XX))
    colnames(out_dat)[1] <- c("Y")
    
    out <- list(list("data_full" = out_dat, "Xmat" = XX, "pscore" = pscore),
                list("influential_vars" = influential_var,"influential_var_indiv" = influential_var_indiv))
  }
  
  return(out)
  
}

#I
f.list.I<-function(var.list, data.x, data.y){
  kk=length(var.list)
  if(kk>1){
    xx=as.matrix(data.x[,as.vector(var.list)]%*%as.vector((3^(0:(kk-1)))))
  }
  else{
    xx=as.matrix(data.x)[,var.list]
  }
  yy=unlist(data.y)
  dat.mat=table(xx,yy)
  n.d=dat.mat[,1]
  n.u=dat.mat[,2]
  nn.d=sum(n.d)
  nn.u=sum(n.u)
  i.score=nn.d*nn.u*sum((n.d/nn.d-n.u/nn.u)^2)/(nn.d+nn.u)
  return(c(var.list, i.score))
}

xdiscrete<-function(x){
  if(length(unique(x))>3){
    temp<-kmeans(x,3,nstart=25)
    res<-temp$cluster-1
  }else{
    res<-x
  }
  return(res)
}

makeDiscrete<-function(X){
  cc <- detectCores()-2
  cl <- makeCluster(cc,type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
  registerDoParallel(cl)
  #cat("Discretizing variables... \n")
  res <-  foreach(xx=1:ncol(X),.packages='Matrix',.export='xdiscrete',.combine='cbind') %dopar% {
    xdiscrete(X[,xx])
  }
  colnames(res)<-colnames(X)
  stopCluster(cl) 
  return(res)
}
#calculate I for marginals all variables
calcI<-function(X,y){
  cc <- detectCores()-2
  cl <- makeCluster(cc,type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
  registerDoParallel(cl)
  res <- foreach(xx=1:ncol(X),.packages='Matrix',.export='f.list.I',.combine=c) %dopar% { 
    out <- f.list.I(1,as.matrix(X[,xx]),as.vector(y))[2]
  }
  stopCluster(cl) 
  names(res)<-colnames(X)
  res<-sort(res,decreasing=TRUE)
  return(res)
}

#calculate I for all 2x all variables
calcI2<-function(X,y,ind){
  cc <- detectCores()-2
  cl <- makeCluster(cc,type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
  registerDoParallel(cl)
  res <- foreach(xx=1:nrow(ind),.packages='Matrix',.export='f.list.I',.combine=c) %dopar% { 
    out <- f.list.I(c(ind[xx,1],ind[xx,2]),as.matrix(X),as.vector(y))[3]
  }
  stopCluster(cl) 
  return(res)
}


## name functions ---------------------------------------------------
##
## this functions is for generating proper file name
##
sim_name <- function(sim_id = NULL, split = NULL, dir = './res'
) {
  
  if (is.null(sim_id) | is.null(split)) {
    stop("please specify all arguments\n")
  }
  
  if (split) {
    sp <- "SP"
  } else {
    sp <- "NS"
  }
  
  
  fn_name_out <- paste(
    dir, "/sim_",
    sim_id, "_",
    sp, ".rds",
    sep = ""
  )
  
  return(fn_name_out)
  
}
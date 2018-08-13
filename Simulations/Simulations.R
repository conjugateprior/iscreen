## Positives ##
#Goal 1: Show differentiating between noisy/influential
#Goal 2: Show finding interactions: a) when marginal present, and b) when marginal not present

## Negatives/Limitations ##
#Goal 3: Discretization of covariates -- can we find the variables?
  # Gaussian
  # Bernoulli
  # Mix
#Goal 4: Full space of covariate interactions -- coverage probas


# **************************************************** #
#                  simulation setup                    #
#                                                      #
#             by: Adeline Lo                           #
#                                                      #
# **************************************************** #
t<-Sys.time()
## Libraries
library(Matrix)
library(Matrix.utils)
require(caret)
library(plyr)
library(dplyr)
library(tidyverse)
library(doParallel)
library(SIS)
#require(randomForest)

## param
n_sim <- 100
#theta <- 1

n_obs <- c(200, 400, 1600)
k_cov <- c(10, 100, 200)
ot_model <- c("LN", "NL", "LNJ", "NLJ", "LNMJ", "NLMJ")
var_type <- c("normal", "bern", "mix")

## full combination
sim_set <- expand.grid(n_obs, k_cov, ot_model, var_type)

# **************************************************** #
#                 read functions                       #
# **************************************************** #

## assume we're in one of three estimator directly
#source("./sim_functions2.R")
source("functions.R")

# **************************************************** #
#                 running simulations                  #
# **************************************************** #

#for (i in 1:nrow(sim_set)) {    ## loop is over different DGP
for (i in 115:126) {    ## loop is over different DGP
  ## local DGP setup
  n_run   <- sim_set[i, 1]      ## number of observations
  k_run   <- sim_set[i, 2]      ## number of covariates
  ot_run  <- sim_set[i, 3]      ## outcome model
  var_run <- sim_set[i, 4]      ## variable type
  
  #res <- rep(NA, n_sim)
  res <- vector("list", n_sim)
  for (g in 1:n_sim) {          ## loop is over iterations (different data given a fixed DGP)
    
    ## generate data
    dat <- sim_dgp(
      nobs = n_run, kcov = k_run,
      outcome = ot_run, #treatment = tr_run,
      var_type = var_run,
      seed = g,
      split = FALSE           ## please change this argument if necessary (if you run sample splitting stuff)
    )
    set.seed(g)
    
    
    ## ###################################### ##
    ## this segment is estimator sepecific    ##
    ## please update the following part       ##
    ## ###################################### ##
    
    ## I screen
    ## Discretize variables if necessary
    XmatD<- makeDiscrete(dat[[1]]$Xmat) 
    
    #Marginals
    iscoreM<-calcI(XmatD,dat[[1]]$data_full[,"Y"])
    iscoreM<-iscoreM[which(iscoreM>0.9)]
    matchedM<-match(dat[[2]]$influential_var_indiv,names(iscoreM))
    M<-dat[[2]]$influential_var_indiv[which(!is.na(matchedM))]
    per.captureM<-length(M)/length(dat[[2]]$influential_var_indiv)
    #Joint
    if(ot_run == "LNJ"|ot_run =="NLJ"|ot_run =="LNMJ"|ot_run =="NLMJ"){
      iscoreJ<-as.data.frame(matrix(NA,ncol=4,nrow=choose(k_run,2)))
      names(iscoreJ)<-c("V1","V2","I","VarSet")
      iscoreJ[,1:2]<-t(combn(k_run,2))
      iscoreJ[,3]<-calcI2(XmatD,dat[[1]]$data_full[,"Y"],iscoreJ)
      iscoreJ<-iscoreJ[which(iscoreJ$I>0.9),]
      iscoreJ<-iscoreJ[order(-iscoreJ$I),]
      iscoreJ$VarSet<-paste("X",iscoreJ$V1,"*X",iscoreJ$V2,sep="")
      matchedJ<-match(dat[[2]]$influential_vars,c(iscoreJ$VarSet,iscoreM))
      J<-dat[[2]]$influential_vars[which(!is.na(matchedJ))]
      per.captureJ<-length(J)/length(dat[[2]]$influential_vars)
    }
    
    ## SIS Screen
    model           <-  SIS(dat[[1]]$Xmat,dat[[1]]$data_full[,"Y"],family="binomial",tune=c("bic"))
    SIS.M           <-  colnames(dat[[1]]$Xmat)[model$ix]
    per.captureSISM <-  length(SIS.M)/length(dat[[2]]$influential_var_indiv)
    
    ## update until this line ############### ##
    ## save relevant info here
    if(ot_run == "LN"|ot_run =="NL"){
    res[[g]]<-list("Captured.M"=M,"Marginal.I"=iscoreM,"Per.captured.M"=per.captureM,
                   "Captured.SIS.M"=SIS.M,"Per.captured.SIS.M"=per.captureSISM)
    }
    else{
    res[[g]]<-list("Captured.M"=M,"Captured.J"=J,"Marginal.I"=iscoreM,"Joint.I"=iscoreJ,
                   "Per.captured.M"=per.captureM, "Per.captured.J"=per.captureJ,
                   "Captured.SIS.M"=SIS.M,"Per.captured.SIS.M"=per.captureSISM)
    }

  }
  
  ##
  ## post processing
  ##
  
  #bias <- mean(res - theta)
  #rmse <- sqrt(mean((res - theta)^2))
  #variance <- var(res)
  
  #% captured sets
  #returned vars/prop. of total variables
  
  
  #sim_out <- c(bias, rmse, variance)
  #names(sim_out) <- c("bias", "rmse", "variance")
  
  ##
  ## save results
  ##
  
  ## ###################################### ##
  ## generate file name                     ##
  ## ###################################### ##
  fn_name <- sim_name(sim_id = i, split = FALSE)

  
  ## ###################################### ##
  
  ## save
  saveRDS(res, file = fn_name)
  cat("Simulation set",i,"of 162 complete.\n")
  
} ## done with all simulations
time<-Sys.time()-t

time
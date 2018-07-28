
# Epanechnikov Kernel function

kernelh <- function(u,h) ifelse(abs(u/h)<1,3/4*(1-(u/h)^2)/h,0)

# Estimate for measurement process: gamma, lambda ----------------------------------------

gammaest <-  function(kerMat,meas_obs_list,covar_list,censor_list,n,p){
  # U gamma wrapper
  part1 <- ugamma1_C(kerMat,covar_list,n,p)
  ugamma <- function(gamma){
    part1 -
      ugamma2_C(gamma,kerMat,meas_obs_list,covar_list,censor_list,n,p)
  }
  #browser()
  gammaestres <- nleqslv(rep(0,p),ugamma)
}

# Estimation of Lambda


lambdaest <-  function(gammahat,kerMat,meas_obs_list,covar_list,censor_list,n,p){
  # U gamma wrapper
  dlambdaest <- dlambda_C(gammahat,kerMat,meas_obs_list,covar_list,censor_list,n,p)
  res <- cbind(meas_obs_list,lapply(dlambdaest, as.vector))
  return(res)
}

# Estimation



#kernel <- function(x,h) exp(-2*abs(x/h))/h
estasy <- function(dataset,kerFun,h,n,p){
  # transform data
  covar_list <- lapply(dataset, function(x) x[["covariates"]])
  meas_obs_list <- lapply(dataset, function(x) x[["meas_times"]])
  obscov_times_list <- lapply(dataset, function(x) x[["obscov_times"]])
  censor_list <- sapply(dataset, function(x) x[["censoring"]])

  # R realization, slow
  #kerMat <- apply(expand.grid(1:n,1:n),1, function(ind) kernelh(outermin_C(dataset[[ind[2]]][[2]],dataset[[ind[1]]][[4]]),h))

  kerMat <- kerMatgen_C(meas_obs_list,obscov_times_list,h)

  #browser()
  gammaest_res <- gammaest(kerMat,meas_obs_list,covar_list,censor_list,n,p)
  # TODO consider add error handling and logging if not converges

  lambdaest_res <- lambdaest(gammaest_res$x,kerMat,meas_obs_list,covar_list,censor_list,n,p)
  Lambdaest_res <- cbind(sort(unlist(lambdaest_res[,1])),cumsum(unlist(lambdaest_res[,2])[order(unlist(lambdaest_res[,1]))]))

  response_list <- sapply(dataset, function(x) x[["Y"]])
  #kerMat <- kerMatgen_C(meas_obs_list,obscov_times_list,nsample ^ (-0.8))

  longest_res <- longest_c(gamma = gammaest_res$x,
                        kerMat=kerMat,
                        meas_times = meas_obs_list,
                        covariates = covar_list,
                        response = response_list,
                        dlambda = NULL,
                        censor = censor_list, n, p)
  #gmu0est <- cbind(unlist(meas_obs_list),unlist(response_list)-do.call(rbind,longest$Xbar_list) %*% longest[[1]])
  #gmu0est <- gmu0est[order(gmu0est[,1]),]
  return(list(gamma=gammaest_res$x,lambda0 = lambdaest_res, Lambda0=Lambdaest_res, thetaest=longest_res[[1]],gmu0est=longest_res[[2]]))
}

################################################################################################################
#kernel <- function(x,h) exp(-2*abs(x/h))/h
estasy_test <- function(dataset,kerFun,h,n,p){
  # transform data
  covar_list <- lapply(dataset, function(x) x[["covariates"]])
  meas_obs_list <- lapply(dataset, function(x) x[["meas_times"]])
  obscov_times_list <- lapply(dataset, function(x) x[["obscov_times"]])
  censor_list <- sapply(dataset, function(x) x[["censoring"]])

  # R realization, slow
  #kerMat <- apply(expand.grid(1:n,1:n),1, function(ind) kernelh(outermin_C(dataset[[ind[2]]][[2]],dataset[[ind[1]]][[4]]),h))

  kerMat <- kerMatgen_C(meas_obs_list,obscov_times_list,h)

  #browser()
  gammaest_res <- gammaest(kerMat,meas_obs_list,covar_list,censor_list,n,p)
  # TODO consider add error handling and logging if not converges

  lambdaest_res <- lambdaest(gammaest_res$x,kerMat,meas_obs_list,covar_list,censor_list,n,p)
  Lambdaest_res <- cbind(sort(unlist(lambdaest_res[,1])),cumsum(unlist(lambdaest_res[,2])[order(unlist(lambdaest_res[,1]))]))


  return(list(gamma=gammaest_res$x,lambda0 = lambdaest_res, Lambda0=Lambdaest_res))
}




# TODO, Do I need to sample censoring time
simdataone <- function(nsample, p, infl, obscov_rate,lambda0_val,mu0,beta,alpha,gamma,censor) {
  simdataoneres <-
    replicate(
      nsample,
      simAsyLongdata(
        obscov_rate = obscov_rate*infl,
        lambda0 = function(x) lambda0_val*infl,
        lmax=NULL,        #lmax = lambda0_val * exp(abs(gamma))*infl,
        nstep=20,
        mu0 = mu0,
        beta = beta,
        alpha = alpha,
        gamma = gamma,
        cen = censor
      ),
      simplify = F
    )
  simdataoneres
}

# clusterEvalQ(cl,{
#   library('asynlong')
#   library("MASS")
#   #library("reda")
#   library("extraDistr")
#   library("nleqslv")
# })

# Main parallel computing
simmain <- function(nrep,nsample, p, infl, obscov_rate,lambda0_val,mu0,beta,alpha,gamma,censor,horder) {
  #cat("Start: ")
  simdatarep <- foreach(i = 1:nrep,.packages=c("asynlong","MASS","extraDistr","nleqslv"),.errorhandling = "remove") %do% {
    simdataone(nsample, p, infl, obscov_rate,lambda0_val,mu0,beta,alpha,gamma,censor)
  }
  #simdatarep <- lapply(1:nrep,function(x) simdataone(nsample, p, infl, obscov_rate,lambda0_val,mu0,beta,alpha,gamma,censor))
  simres <- foreach(simdata = simdatarep,.combine="rbind",.packages=c("asynlong","MASS","extraDistr","nleqslv"),.errorhandling = "remove") %dopar% {

    simoneres <- estasy_pur(simdata, NULL, nsample ^ (-horder), nsample, p)
    CI_res <- simoneres$CI_theta
    Cp_ind <- (CI_res[,1] < c(gamma,beta,alpha)) * (CI_res[,2] > c(gamma,beta,alpha))
    CI_res_pur <- simoneres$CI_pur
    Cp_ind_pur <- (CI_res_pur[,1] < c(gamma,beta,alpha)) * (CI_res_pur[,2] > c(gamma,beta,alpha))
    c(simoneres[[1]],as.vector(Cp_ind),as.vector(t(CI_res)),as.vector(Cp_ind_pur),as.vector(t(CI_res_pur)))
  }
  #print(simres)
  #cat("End: ")
  return(simres)
}


simrun <- function(nrep,nsample,p,infl=2,obscov_rate,lambda0_val,mu0,beta0,alpha0,gamma0,censor=1,horder) {
  cat("PID: ", Sys.getpid(), " Start: ",format(Sys.time(), "%a %b %d %X %Y"))
  simres <- simmain(
    nrep = nrep,
    nsample = nsample,
    p = p,
    infl = infl,
    obscov_rate = obscov_rate,
    lambda0_val = lambda0_val,
    mu0 = mu0,
    beta = beta0,
    alpha = alpha0,
    gamma = gamma0,
    censor = censor,
    horder = horder
  )
  cat("End: \n")
  poiest <- simres[,1:length(c(gamma0,beta0,alpha0))]
  bias <- colMeans(poiest)-c(gamma0,beta0,alpha0)
  variance <- diag(var(poiest))
  stdev <- sqrt(variance)
  rmse <- sqrt(colMeans(t((t(poiest)-c(gamma0,beta0,alpha0))^2)))
  Cp <- colMeans(simres[,-c(1:length(c(gamma0,beta0,alpha0)))])
  data.frame(bias,variance,stdev,rmse,Cp=Cp[1:length(c(gamma0,beta0,alpha0))])
}



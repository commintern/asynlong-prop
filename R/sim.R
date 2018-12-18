
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


# Main parallel computing
simmain <- function(nrep,nsample, p, infl, obscov_rate,lambda0_val,mu0,beta,alpha,gamma,censor,horder) {
  #cat("Start: ")
  simdatarep <- parLapply(cl,1:nrep,function(x) simdataone(nsample, p, infl, obscov_rate,lambda0_val,mu0,beta,alpha,gamma,censor))
  #simdatarep <- lapply(1:nrep,function(x) simdataone(nsample, p, infl, obscov_rate,lambda0_val,mu0,beta,alpha,gamma,censor))
  simres <- foreach(simdata = simdatarep,.combine="rbind") %dopar% {

    simoneres <- estasy(simdata, NULL, nsample ^ (-horder), nsample, p)
    unlist(simoneres[c(1,4)])
  }
  #cat("End: ")
  return(simres)
}


simrun <- function(nrep,nsample,p,infl=2,obscov_rate,lambda0_val,mu0,beta0,alpha0,gamma0,censor=1,horder) {
  cat("Start: horder is",horder)
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

  bias <- colMeans(simres)-c(gamma0,beta0,alpha0)
  variance <- diag(var(simres))
  stdev <- sqrt(variance)
  rmse <- sqrt(colMeans(t((t(simres)-c(gamma0,beta0,alpha0))^2)))
  data.frame(bias,variance,stdev,rmse)
}



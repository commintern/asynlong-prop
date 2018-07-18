
# TODO, Do I need to sample censoring time
simdataone <- function(nsample, p, infl, obscov_rate,lambda0_val,mu0,beta,alpha,gamma,censor) {
  simdatasam <-
    replicate(
      nsample,
      simAsyLongdata(
        obscov_rate = obscov_rate*infl,
        lambda0 = function(x) lambda0_val*infl,
        lmax = lambda0_val * exp(abs(gamma))*infl,
        mu0 = mu0,
        beta = beta,
        alpha = alpha,
        gamma = gamma,
        cen = censor
      ),
      simplify = F
    )
  simdatasam
}


simmain <- function(nrep,nsample, p, infl, obscov_rate,lambda0_val,mu0,beta,alpha,gamma,cenor,horder) {
  simdatarep <- replicate(nrep,simdataone(nsample, p, infl, obscov_rate,lambda0_val,mu0,beta,alpha,gamma,cenor))
  simres <- foreach(simdata = simdatarep,.combine="rbind") %do% {
    simoneres <- estasy(simdata, NULL, nsample ^ (-horder), nsample, p)
    unlist(simoneres[c(1,4)])
  }
  simres
}


# simone <- function(nsample, p, infl, obscov_rate,lambda0,mu0,beta,alpha,gamma,cen,horder) {
#   simdatasam <-
#     replicate(
#       nsample,
#       simAsyLongdata(
#         obscov_rate = 6*infl,
#         lambda0 = function(x) 2*infl,
#         lmax = 2 * exp(1)*infl,
#         mu0 = function(x) exp(sin(2 * pi * x)),
#         beta = 2,
#         alpha = 1,
#         gamma = 1,
#         cen = 1
#       ),
#       simplify = F
#     )
#   simoneres <- estasy(simdatasam, NULL, nsample ^ (-0.8), nsample, p)
#   simoneres
# }

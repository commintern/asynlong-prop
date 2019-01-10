if(!require("doParallel")) install.packages("doParallel",repos = "http://cran.us.r-project.org")
library("doParallel")
 library('asynlong')
#stopCluster(cl)

cl <- makeCluster(4,outfile="debug.txt")
registerDoParallel(cl)
clusterEvalQ(cl,{
  library('asynlong')
  library("MASS")
  #library("reda")
  library("extraDistr")
  library("nleqslv")
})
#clusterCall(cl, library, package=c('asynlong'))


simcol <- function(horderlist,nrep,nsample,p,infl,obscov_rate,lambda0_val,mu0,beta0,alpha0,gamma0,censor) {
  res <- sapply(horderlist, function(x) simrun(nrep = nrep, nsample = nsample,p = p,infl = infl, obscov_rate = obscov_rate,lambda0_val = lambda0_val,
                                               mu0 = mu0, beta0 = beta0, alpha0 = alpha0, gamma0 = gamma0, censor = censor, horder = x),simplify = F)

  res <- do.call(rbind,mapply(function(A,x) cbind(A,x),res,horderlist,SIMPLIFY=F))
  l <- length(horderlist)

  cbind(nsample,p,infl,para=rep(c("gamma","beta","alpha"),l),tr=rep(c(gamma0,beta0,alpha0),l),res)
}

# simcol(nrep = 1000, nsample = 100,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = c(0.6,0.7,0.8,0.9))
#
# simcol(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = c(0.6,0.7,0.8,0.9))
#
# simcol(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = c(0.6,0.7,0.8,0.9))
#
# simcol(nrep = 1000, nsample = 2000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = c(0.6,0.7,0.8,0.9))
#
# simcol(nrep = 1000, nsample = 3000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = c(0.6,0.7,0.8,0.9))


#######################################################################################



simcol(nrep = 500, nsample = 100,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 6,
       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = c(6:8)/10)

simcol(nrep = 1000, nsample = 500,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = c(2:8)/10)

simcol(nrep = 1000, nsample = 1000,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = c(2:8)/10)

simcol(nrep = 1000, nsample = 100,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 1, censor = 1, horder = c(2:8)/10)

simcol(nrep = 1000, nsample = 500,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 1, censor = 1, horder = c(2:8)/10)

simcol(nrep = 1000, nsample = 1000,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 1, censor = 1, horder = c(2:8)/10)

simcol(nrep = 1000, nsample = 100,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 0, censor = 1, horder = c(2:8)/10)

simcol(nrep = 1000, nsample = 500,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 0, censor = 1, horder = c(2:8)/10)

simcol(nrep = 1000, nsample = 1000,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 0, censor = 1, horder = c(2:8)/10)

simcol(nrep = 1000, nsample = 100,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 0, censor = 1, horder = c(2:8)/10)

simcol(nrep = 1000, nsample = 500,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 0, censor = 1, horder = c(2:8)/10)

simcol(nrep = 1000, nsample = 1000,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 0, censor = 1, horder = c(2:8)/10)

#simcol(nrep = 1000, nsample = 1600,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
#       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = c(2:8)/10)

#simcol(nrep = 1000, nsample = 3000,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
#       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = c(2:5)/10)


# simrun(nrep = 1000, nsample = 100,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.9)
#
# simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.9)
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#   mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.9)
#
#
# simrun(nrep = 1000, nsample = 2000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.9)
#
# simrun(nrep = 1000, nsample = 5000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.9)
#
#
#
# simrun(nrep = 1000, nsample = 100,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.8)
#
# simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.8)
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.8)
#
#
# simrun(nrep = 1000, nsample = 2000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.8)
#
# simrun(nrep = 1000, nsample = 5000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.8)
#
#
# ##############################################################################################################
#
#
# simrun(nrep = 1000, nsample = 100,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.9)
#
# simrun(nrep = 1000, nsample = 500,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.9)
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.9)
#
#
# simrun(nrep = 1000, nsample = 2000,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.9)
#
# simrun(nrep = 1000, nsample = 5000,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.9)
#
#
#
# simrun(nrep = 1000, nsample = 100,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.8)
#
# simrun(nrep = 1000, nsample = 500,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.8)
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.8)
#
#
# simrun(nrep = 1000, nsample = 2000,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.8)
#
# simrun(nrep = 1000, nsample = 5000,p = 1,infl = 2, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.8)
#

# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#   mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 1, censor = 1, horder = 0.6)
#
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 1, censor = 1, horder = 0.7)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 1, censor = 1, horder = 0.6)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 1, censor = 1, horder = 0.7)
#
# #=-=====================================================================================================
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.6)
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.7)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.6)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.7)
#
#
# #==================================================================================================
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 0, censor = 1, horder = 0.6)
#
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 0, censor = 1, horder = 0.7)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 0, censor = 1, horder = 0.6)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 0, gamma0 = 0, censor = 1, horder = 0.7)
#
#
# #=-=====================================================================================================
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 0, censor = 1, horder = 0.6)
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 0, censor = 1, horder = 0.7)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 0, censor = 1, horder = 0.6)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sin(2 * pi * x)), beta0 = -2, alpha0 = 1, gamma0 = 0, censor = 1, horder = 0.7)
#
# ##################################################################################################################
# ######################################################################################################################
#
# #==================================================================================================
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 0, gamma0 = 1, censor = 1, horder = 0.6)
#
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 0, gamma0 = 1, censor = 1, horder = 0.7)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 0, gamma0 = 1, censor = 1, horder = 0.6)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 0, gamma0 = 1, censor = 1, horder = 0.7)
#
# #=-=====================================================================================================
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.6)
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.7)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.6)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 1, gamma0 = 1, censor = 1, horder = 0.7)
#
#
# #==================================================================================================
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 0, gamma0 = 0, censor = 1, horder = 0.6)
#
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 0, gamma0 = 0, censor = 1, horder = 0.7)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 0, gamma0 = 0, censor = 1, horder = 0.6)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 0, gamma0 = 0, censor = 1, horder = 0.7)
#
#
# #=-=====================================================================================================
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 1, gamma0 = 0, censor = 1, horder = 0.6)
#
# simrun(nrep = 1000, nsample = 1000,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
#        mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 1, gamma0 = 0, censor = 1, horder = 0.7)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 1, gamma0 = 0, censor = 1, horder = 0.6)
#
# #simrun(nrep = 1000, nsample = 500,p = 1,infl = 1, obscov_rate = 6,lambda0_val = 2,
# #       mu0 = function(x) exp(sqrt(x)), beta0 = -2, alpha0 = 1, gamma0 = 0, censor = 1, horder = 0.7)

stopCluster(cl)

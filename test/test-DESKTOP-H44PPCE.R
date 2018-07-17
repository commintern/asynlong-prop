library("nleqslv")
library("MASS")
library("reda")
library("extraDistr")


# test of data geneation--------------------------------------------------------------------

# univariate z
simAsyLongdata(obscov_rate=5, lambda0=function(x) 5, lmax=5*exp(1.5), mu0=function(x) exp(sin(2*pi*x)), beta=1.5, alpha=1, gamma=1,cen=1)

simdatasam <- replicate(500,simAsyLongdata(obscov_rate=5, lambda0=function(x) 5, lmax=5*exp(1.5), mu0=function(x) exp(sin(2*pi*x)), beta=1.5, alpha=1, gamma=1,cen=1)
                        ,simplify = F)

# multivairate z

simAsyLongdata(obscov_rate=5, lambda0=function(x) 5, lmax=5*exp(1.5), mu0=function(x) exp(sin(2*pi*x)), beta=c(0.5,1.5), alpha=1, gamma=c(0.5,1),cen=1)

simdatasam <- replicate(10,simAsyLongdata(obscov_rate=5, lambda0=function(x) 5, lmax=5*exp(1.5), mu0=function(x) exp(sin(2*pi*x)), beta=c(0.5,1.5), alpha=1, gamma=c(0.5,1),cen=1)
                        ,simplify = F)
#simdatasam <- replicate(10,simAsyLongdata(5, 5, 0.5, 1.5, 2, 1))

#simdatasam <- replicate(10,simAsyLongdata(5, 5, 0.5, 1.5, 2, 1),simplify=F)

testcount <- function(n){
  t1 <- sort(runif(1 + rpois(1,n)))
  t2 <- sort(runif(1 + rpois(1,n)))
  stepfun(t1,0:length(t1),right=TRUE)(t2)-as.vector(countprofun_C(t1,t2))
}


# Test of obsevration process estimation ----------------------------------------------------

source('datagen.R')
source('estimation.R')
Rcpp::sourceCpp('Cpp/util.cpp')
Rcpp::sourceCpp('Cpp/recest.cpp')
Rcpp::sourceCpp('Cpp/longest.cpp')

library("nleqslv")
library("MASS")
library("reda")
library("extraDistr")

set.seed(8102)
nsample <- 500

simdatasam <-
  replicate(
    nsample,
    simAsyLongdata(
      obscov_rate = 5,
      lambda0 = function(x) 5,
      lmax = 5 * exp(2),
      mu0 = function(x) exp(sin(2 * pi * x)),
      beta = 1.5,
      alpha = 1,
      gamma = 2,
      cen = 1
    ),
    simplify = F
  )
testres <- estasy(simdatasam, kernelh, nsample ^ (-0.5), nsample, 1)
testres[[1]]
plot(testres[[3]], typ = "l")
abline(0, 5, col = "red")

# Test of longitudinal data estimation

covar_list <- lapply(simdatasam, function(x) x[["covariates"]])
meas_obs_list <- lapply(simdatasam, function(x) x[["meas_times"]])
obscov_times_list <- lapply(simdatasam, function(x) x[["obscov_times"]])
censor_list <- sapply(simdatasam, function(x) x[["censoring"]])
dlambda_list <-  testres[[2]][,2]
response_list <- sapply(simdatasam, function(x) x[["Y"]])
kerMat <- kerMatgen_C(meas_obs_list,obscov_times_list,nsample ^ (-0.9))

testlong <- longest_c(gamma = 2,
          kerMat=kerMat,
          meas_times = meas_obs_list,
          covariates = covar_list,
          response = response_list,
          dlambda = dlambda_list,
          censor = censor_list, nsample, 1) 

testlong[[1]]

testlong <- longest_test_c(gamma = 2,
                      kerMat=kerMat,
                      meas_times = meas_obs_list,
                      covariates = covar_list,
                      response = response_list,
                      dlambda = dlambda_list,
                      censor = censor_list, nsample, 1) 

testlong[[1]]

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
nsample <- 100


infl <- 2
simdatasam <-
  replicate(
    nsample,
    simAsyLongdata(
      obscov_rate = 6*infl,
      lambda0 = function(x) 2*infl,
      lmax = 2 * exp(1)*infl,
      mu0 = function(x) exp(sin(2 * pi * x)),
      beta = 2,
      alpha = 1,
      gamma = 1,
      cen = 1
    ),
    simplify = F
  )
testres <- estasy_test(simdatasam, kernelh, nsample ^ (-0.8), nsample, 1)
testres[[1]]
plot(testres[[3]], typ = "l")
abline(0, 2*infl, col = "red")

# Test of longitudinal data estimation

covar_list <- lapply(simdatasam, function(x) x[["covariates"]])
meas_obs_list <- lapply(simdatasam, function(x) x[["meas_times"]])
obscov_times_list <- lapply(simdatasam, function(x) x[["obscov_times"]])
censor_list <- sapply(simdatasam, function(x) x[["censoring"]])
dlambda_list <-  testres[[2]][,2]
response_list <- sapply(simdatasam, function(x) x[["Y"]])
kerMat <- kerMatgen_C(meas_obs_list,obscov_times_list,nsample ^ (-0.8))

testlong <- longest_c(gamma = testres[[1]],
          kerMat=kerMat,
          meas_times = meas_obs_list,
          covariates = covar_list,
          response = response_list,
          dlambda = dlambda_list,
          censor = censor_list, nsample, 1) 

testlong[[1]]

testlong <- longest_test_c(gamma = testres[[1]],
                      kerMat=kerMat,
                      meas_times = meas_obs_list,
                      covariates = covar_list,
                      response = response_list,
                      dlambda = dlambda_list,
                      censor = censor_list, nsample, 1) 

testlong[[1]]




infl <- 2
p <- 2
simdatasam <-
  replicate(
    nsample,
    simAsyLongdata(
      obscov_rate = 6*infl,
      lambda0 = function(x) 2*infl,
      lmax = 2 * exp(1)*infl,
      mu0 = function(x) exp(sin(2 * pi * x)),
      beta = c(2,4),
      alpha = 0,
      gamma = c(1,1)/2,
      cen = 1
    ),
    simplify = F
  )
testres <- estasy_test(simdatasam, kernelh, nsample ^ (-0.8), nsample, p)
testres[[1]]
plot(testres[[3]], typ = "l")
abline(0, 2*infl, col = "red")

# Test of longitudinal data estimation

covar_list <- lapply(simdatasam, function(x) x[["covariates"]])
meas_obs_list <- lapply(simdatasam, function(x) x[["meas_times"]])
obscov_times_list <- lapply(simdatasam, function(x) x[["obscov_times"]])
censor_list <- sapply(simdatasam, function(x) x[["censoring"]])
dlambda_list <-  testres[[2]][,2]
response_list <- sapply(simdatasam, function(x) x[["Y"]])
kerMat <- kerMatgen_C(meas_obs_list,obscov_times_list,nsample ^ (-0.8))

testlong <- longest_c(gamma = testres[[1]],
                      kerMat=kerMat,
                      meas_times = meas_obs_list,
                      covariates = covar_list,
                      response = response_list,
                      dlambda = dlambda_list,
                      censor = censor_list, nsample, p) 

testlong[[1]]


###################################################################################################

simone <- function(nsample, infl, p) {
  simdatasam <-
    replicate(
      nsample,
      simAsyLongdata(
        obscov_rate = 6*infl,
        lambda0 = function(x) 2*infl,
        lmax = 2 * exp(1)*infl,
        mu0 = function(x) exp(sin(2 * pi * x)),
        beta = 2,
        alpha = 1,
        gamma = 1,
        cen = 1
      ),
      simplify = F
    )
  simoneres <- estasy(simdatasam, NULL, nsample ^ (-0.8), nsample, p)
  simoneres
}

tt <- simone(100,1,1)

system.time( simone(500,1,1))

library(doParallel)

stopCluster(cl)

cl <- cl <- makeCluster(6, type = "PSOCK") 
registerDoParallel(cl)
clusterExport(cl,"simone")
clusterEvalQ(cl,{
  source('datagen.R')
  source('estimation.R')
  Rcpp::sourceCpp('Cpp/util.cpp')
  Rcpp::sourceCpp('Cpp/recest.cpp')
  Rcpp::sourceCpp('Cpp/longest.cpp')
  
  library("nleqslv")
  library("MASS")
  library("reda")
  library("extraDistr")
},env=.GlobalEnv)

clusterEvalQ(cl,{
  source('datagen.R')
  source('estimation.R')
})

clusterCall(cl, Rcpp::sourceCpp, file='Cpp/util.cpp',env=.GlobalEnv)
clusterCall(cl, epanker_C,0.5,1)


parSapply(cl = NULL, X, FUN, ..., simplify = TRUE,
          USE.NAMES = TRUE, chunk.size = NULL)

foreach(i=1:6) %dopar% {
  1+1
}



simres <- foreach(i=1:1000) %do% {
  simresone <- simone(100,1,1)
  simresone[c(1,4)]
}
simres <- t(sapply(simres, function(x) c(x[[1]],x[[2]])))

colMeans(simres)

simres <- foreach(i=1:1000) %do% {
  simresone <- simone(200,1,1)
  simresone[c(1,4)]
}


simres <- t(sapply(simres, function(x) c(x[[1]],x[[2]])))

colMeans(simres)

simres <- foreach(i=1:100) %do% {
  simresone <- simone(500,1,1)
  simresone[c(1,4)]
}
simres <- t(sapply(simres, function(x) c(x[[1]],x[[2]])))

colMeans(simres)
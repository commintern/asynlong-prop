library(doParallel)

#stopCluster(cl)

cl <- makeCluster(24,outfile="debug.txt")
registerDoParallel(cl)
clusterEvalQ(cl,{
  library('asynlong')
  library("MASS")
  #library("reda")
  library("extraDistr")
  library("nleqslv")
})
#clusterCall(cl, library, package=c('asynlong'))

simrun(
  nrep = 100,
  nsample = 100,
  p = 1,
  infl = 2,
  obscov_rate = 6,
  lambda0_val = 2,
  mu0 = function(x)
    exp(sin(2 * pi * x)),
  beta0 = 2,
  alpha0 = 1,
  gamma0 = 1,
  censor = 1,
  horder = 0.6
)

stopCluster(cl)

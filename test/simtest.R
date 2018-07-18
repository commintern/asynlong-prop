library(doParallel)

stopCluster(cl)

cl <- makeCluster(4, type = "PSOCK")
registerDoParallel(cl)
clusterEvalQ(cl,{
  library('asynlong')
  library("MASS")
  #library("reda")
  library("extraDistr")
  library("nleqslv")
})
#clusterCall(cl, library, package=c('asynlong'))



simmain(
  nrep = 2,
  nsample = 100,
  p = 1,
  infl = 2,
  obscov_rate = 6,
  lambda0_val = 2,
  mu0 = function(x)
    exp(sin(2 * pi * x)),
  beta = 2,
  alpha = 1,
  gamma = 1,
  cenor = 1,
  horder = -0.8
)



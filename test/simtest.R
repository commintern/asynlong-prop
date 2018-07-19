library(doParallel)

stopCluster(cl)

cl <- makeCluster(6, type = "PSOCK")
registerDoParallel(cl)
clusterEvalQ(cl,{
  library('asynlong')
  library("MASS")
  #library("reda")
  library("extraDistr")
  library("nleqslv")
})
#clusterCall(cl, library, package=c('asynlong'))



simres <- simmain(
  nrep = 500,
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
  horder = 0.6
)

colMeans(simres)
var(simres)
sqrt(diag(var(simres)))
sqrt(colMeans(t((t(simres)-c(1,2,1))^2)))




simres <- simmain(
  nrep = 100,
  nsample = 300,
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
  horder = 0.9
)

colMeans(simres)
var(simres)



simres <- simmain(
  nrep = 100,
  nsample = 100,
  p = 1,
  infl = 2,
  obscov_rate = 6,
  lambda0_val = 2,
  mu0 = function(x) sqrt(x),
  beta = 2,
  alpha = 1,
  gamma = 1,
  cenor = 1,
  horder = 0.4
)
colMeans(t((t(simres)-c(1,2,1))^2))



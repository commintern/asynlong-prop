library(doParallel)

stopCluster(cl)

cl <- makeCluster(6)
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
  nrep = 100,
  nsample = 100,
  p = 1,
  infl = 2,
  obscov_rate = 6,
  lambda0_val = 4,
  mu0 = function(x)
    exp(sin(2 * pi * x)),
  beta = 2,
  alpha = 1,
  gamma = 1,
  censor = 1,
  horder = 0.6
)

colMeans(simres)
var(simres)
sqrt(diag(var(simres)))
sqrt(colMeans(t((t(simres)-c(1,2,1))^2)))




simres <- simmain(
  nrep = 100,
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
  mu0 = function(x)
    sqrt(x),
  beta = 2,
  alpha = 1,
  gamma = 1,
  censor = 1,
  horder = 0.4
)
colMeans(t((t(simres)-c(1,2,1))^2))


library(doParallel)

stopCluster(cl)

cl <- makeCluster(6)
registerDoParallel(cl)
clusterEvalQ(cl,{
  library('asynlong')
  library("MASS")
  #library("reda")
  library("extraDistr")
  library("nleqslv")
})

simrun(
  nrep = 1000,
  nsample = 900,
  p = 1,
  infl = 1,
  obscov_rate = 6,
  lambda0_val = 3,
  mu0 = function(x)
    exp(2),
  beta0 = 2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.6
)

simrun(
  nrep = 1000,
  nsample = 400,
  p = 1,
  infl = 1,
  obscov_rate = 6,
  lambda0_val = 3,
  mu0 = function(x)
    exp(2),
  beta0 = 2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.6
)

simrun(
  nrep = 10000,
  nsample = 100,
  p = 1,
  infl = 1,
  obscov_rate = 8,
  lambda0_val = 4,
  mu0 = function(x)
    exp(2),
  beta0 = -2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.5
)



simtemp <- function(nrep,nsample,p,infl=2,obscov_rate,lambda0_val,mu0,beta0,alpha0,gamma0,censor=1,horder) {
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
  simres
}


r1p <- simtemp(
  nrep = 100,
  nsample = 100,
  p = 1,
  infl = 1,
  obscov_rate = 8,
  lambda0_val = 4,
  mu0 = function(x)
    exp(2),
  beta0 = -2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.66
)

r1p2 <- simtemp(
  nrep = 400,
  nsample = 100,
  p = 1,
  infl = 1,
  obscov_rate = 8,
  lambda0_val = 4,
  mu0 = function(x)
    exp(2),
  beta0 = -2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.66
)

r1p3 <- simtemp(
  nrep = 10,
  nsample = 50,
  p = 1,
  infl = 1,
  obscov_rate = 8,
  lambda0_val = 4,
  mu0 = function(x)
    exp(2),
  beta0 = -2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.5
)

r1 <- simtemp(
  nrep = 10,
  nsample = 100,
  p = 1,
  infl = 1,
  obscov_rate = 8,
  lambda0_val = 4,
  mu0 = function(x)
    exp(2),
  beta0 = -2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.66
)
colMeans(r1)

 r11 <- simtemp(
  nrep = 1000,
  nsample = 100,
  p = 1,
  infl = 1,
  obscov_rate = 8,
  lambda0_val = 4,
  mu0 = function(x)
    exp(2),
  beta0 = -2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.6
)
colMeans(r11)

r111 <- simtemp(
  nrep = 1000,
  nsample = 100,
  p = 1,
  infl = 1,
  obscov_rate = 10,
  lambda0_val = 5,
  mu0 = function(x)
    exp(2),
  beta0 = -2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.66
)
colMeans(r111)

r112 <- simtemp(
  nrep = 2000,
  nsample = 200,
  p = 1,
  infl = 1,
  obscov_rate = 10,
  lambda0_val = 5,
  mu0 = function(x)
    exp(2),
  beta0 = -2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.66
)
colMeans(r112)



r12 <- simtemp(
  nrep = 1000,
  nsample = 200,
  p = 1,
  infl = 1,
  obscov_rate = 8,
  lambda0_val = 4,
  mu0 = function(x)
    exp(2),
  beta0 = -2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.6
)
colMeans(r12)

r2 <- simtemp(
  nrep = 500,
  nsample = 300,
  p = 1,
  infl = 1,
  obscov_rate = 8,
  lambda0_val = 4,
  mu0 = function(x)
    exp(2),
  beta0 = -2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.6
)
colcolMeans(r2)

r3 <- simtemp(
  nrep = 1000,
  nsample = 900,
  p = 1,
  infl = 1,
  obscov_rate = 8,
  lambda0_val = 4,
  mu0 = function(x)
    exp(2),
  beta0 = -2,
  alpha0 = 1,
  gamma0 = 1.5,
  censor = 1,
  horder = 0.66
)
 colMeans(r3)


 r31 <- simtemp(
   nrep = 1000,
   nsample = 900,
   p = 1,
   infl = 1,
   obscov_rate = 8,
   lambda0_val = 4,
   mu0 = function(x)
     exp(2),
   beta0 = -2,
   alpha0 = 1,
   gamma0 = 1.5,
   censor = 1,
   horder = 0.7
 )
 colMeans(r31)

 r4 <- simtemp(
   nrep = 50,
   nsample = 1500,
   p = 1,
   infl = 1,
   obscov_rate = 6,
   lambda0_val = 3,
   mu0 = function(x)
     exp(2),
   beta0 = -2,
   alpha0 = 1,
   gamma0 = 1.5,
   censor = 1,
   horder = 0.6
 )
 colMeans(r4)

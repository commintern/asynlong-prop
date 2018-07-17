library("MASS")
tt <- function(n) {
  ts <- seq(0,1,length.out=n)
  Sigtest <- exp(-abs(outer(ts,ts,"-")))
  test <- mvrnorm(1,rep(0,n),Sigtest)
  diff(ts)[1]*sum(exp(test))
}

tr <- replicate(100,tt(1000))
mean(tr)
tr <- tt()



kerMat[[1]] %*% (t(simdatasam[[1]]$covariates) * exp(colSums(simdatasam[[1]]$covariates)))+
  kerMat[[2]] %*% (t(simdatasam[[2]]$covariates) * exp(colSums(simdatasam[[2]]$covariates)))

kerMat[[1]] %*% ( exp(colSums(simdatasam[[1]]$covariates)))+
  kerMat[[2]] %*% (  exp(colSums(simdatasam[[2]]$covariates)))

simdatasam <- replicate(2,simAsyLongdata(5, 5, 0.5, 1.5, 2, 1),simplify=F)
kerMat <- apply(expand.grid(1:2,1:2),1, function(ind) kernel(outer(simdatasam[[ind[2]]][[2]],simdatasam[[ind[1]]][[4]],"-"),1))

covarlist <- lapply(simdatasam, function(x) x[["covariates"]])
ugamma1_C(kerMat,covarlist,2L,1)

library("MASS")
simdatasam <- replicate(2,simAsyLongdata(5, 5, 0.5, 1.5, 2, 1),simplify=F)
kernel <- function(x,h) exp(-2*abs(x/h))/h
kerMat <- apply(expand.grid(1:2,1:2),1, function(ind) kernel(outer(simdatasam[[ind[2]]][[2]],simdatasam[[ind[1]]][[4]],"-"),1))

covarlist <- lapply(simdatasam, function(x) x[["covariates"]])
timelist <- lapply(simdatasam, function(x) x[["meas_times"]])
censorlist <- sapply(simdatasam, function(x) x[["censoring"]])
ugamma1_C(kerMat,covarlist,2,1)

zbar(1,kerMat,timelist,covarlist,censorlist,2,1)

simdatasam <- replicate(2,simAsyLongdata(5, 5, 0.5, c(1,1.5), 2, 1),simplify=F)
kernel <- function(x,h) exp(-2*abs(x/h))/h
kerMat <- apply(expand.grid(1:2,1:2),1, function(ind) kernel(outer(simdatasam[[ind[2]]][[2]],simdatasam[[ind[1]]][[4]],"-"),1))

covarlist <- lapply(simdatasam, function(x) x[["covariates"]])
timelist <- lapply(simdatasam, function(x) x[["meas_times"]])
censorlist <- sapply(simdatasam, function(x) x[["censoring"]])
ugamma1_C(kerMat,covarlist,2,2)

zbar(c(1,1),kerMat,timelist,covarlist,censorlist,2,2)

ugamma1_C(kerMat,covarlist,2,2)-ugamma2_C(c(1,1),kerMat,timelist,covarlist,censorlist,2,2)

microbenchmark(R=rowSums(do.call(cbind,lapply(simdatasam, function(x) ugam1ind(x,kernel)))),
               Rcpp=ugamma1_C(kerMat,covarlist,10L,2),times=10L)



simdatasam <- replicate(100,simAsyLongdata(5, 5, 0.5, 1.5, 2, 1),simplify=F)
kernel <- function(x,h) exp(-2*abs(x/h))/h
kerMat <- apply(expand.grid(1:100,1:100),1, function(ind) kernel(outer(simdatasam[[ind[2]]][[2]],simdatasam[[ind[1]]][[4]],"-"),1))

covarlist <- lapply(simdatasam, function(x) x[["covariates"]])
timelist <- lapply(simdatasam, function(x) x[["meas_times"]])
censorlist <- sapply(simdatasam, function(x) x[["censoring"]])
ugamma <- function(gamma) ugamma1_C(kerMat,covarlist,100,1)-ugamma2_C(gamma,kerMat,timelist,covarlist,censorlist,100,1)

S0_C(c(1),kerMat,timelist,covarlist,censorlist,100,1)
library("nleqslv")
microbenchmark(nleqslv(c(0),ugamma))


###################
# Test reda package simulation homogenouse poisson process

as.vector(simEvent(rho = 1, origin = 0, endTime = 5))

simhomoPoipro <- function(rate,endTime){
  num <- rpois(1,rate*endTime)
  return(endTime*runif(num))
}


microbenchmark::microbenchmark(as.vector(simEvent(rho = 1, origin = 0, endTime = 5)),simhomoPoipro(1,5))





data <- simdatasam[[1]]
foreach(data=simdatasam,.combine='+') %do% {
foreach(time = data$meas_times, .combine='+') %do% {
  weightsum <- sum(kernelh(time-data$obscov_times,2))
  if(weightsum==0){
    res <-0
  } else {
    res <- data$covariates %*% kernelh(time-data$obscov_times,2)/weightsum
  }
  as.vector(res)
}
}
data$covariates %*% kernelh(t-data$obscov_times,2)



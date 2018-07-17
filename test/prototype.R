# Assume intensity lambda(t)=exp(beta*z(t)), z(t)= 1/(1+exp(x(t)))
library("MASS")
# Assume beta=2
recgaupro <- function(beta) {
  lambdamax <- exp(beta)
  
  tau <- 5
  npoint <- 0
  while(npoint==0){
    npoint <- rpois(1,lambdamax*tau)
  }
  
  selength <-0
  while(selength ==0){
    tall <- sort(runif(npoint)*tau)
    Sigmat <- 2^(-abs(outer(tall,tall,"-")))
    xori <- mvrnorm(n=1,mu=rep(0,length(tall)),Sigma=Sigmat)
    z <- 1/(1+exp(xori))
    lambdat <- exp(beta*z)
    selidx <- runif(length(tall))<lambdat/lambdamax
    selength <- length(selidx)
    tgen <- tall[selidx]
    zgen <- z[selidx]
  }
  
  return(cbind(tgen,zgen))
}

datasetgen <- function(beta,id){
  oridata <- recgaupro(beta)
  temp <- cbind(id,oridata,1,0)
  colnames(temp) <- c("id","time","cov","event","death")
  return(temp)
}



simdatares <- data.frame(do.call(rbind,sapply(1:100,function(id) datasetgen(0.5,id))))
#reReg(reSurv(t,event,death,id)~cov,data=simdatares,method == "cox.LWYY")
fit <-coxph(Surv(t,event)~cov,simdatares)
plot(survfit(fit))

rateReg(Survr(id,time,event)~cov,data=simdatares)


simdatares <- sapply(1:100,function(id) datasetgen(1,id))

kh <- function(u,h) ifelse(abs(u/h)<1,3/4*(1-(u/h)^2)/h,0)

lamesttest <- function(t,beta,simdatares,n,h){
  s0 <-0
  for(i in 1:n){
    tl <- length(simdatares[[i]][,"time"])
    weightsum <- sum(kh(t-simdatares[[i]][-tl,"time"],h))
    if (weightsum!=0){
      s0<- s0+sum(kh(t-simdatares[[i]][-tl,"time"],h)*exp(beta*simdatares[[i]][-tl,"cov"]))/sum(kh(t-simdatares[[i]][-tl,"time"],h))
      #s0<- s0+sum(kh(t-simdatares[[i]][-tl,"time"],1)*exp(beta*simdatares[[i]][-tl,"cov"]))
    }
    
    
  }
  return(1/s0)
}

lamesttest <- function(t,beta,simdatares,n,h){
  s0 <-0
  temp <-0
  for(i in 1:n){
    tl <- length(simdatares[[i]][,"time"])
    if(t %in% simdatares[[i]][-tl,"time"]){
      weightsum <- sum(kh(t-simdatares[[i]][-tl,"time"],h))
      temp <- temp + weightsum
    }
    
    
      s0<- s0+sum(kh(t-simdatares[[i]][-tl,"time"],h)*exp(beta*simdatares[[i]][-tl,"cov"]))
      
    
    
    #s0<- s0+sum(kh(t-simdatares[[i]][-tl,"time"],1)*exp(beta*simdatares[[i]][-tl,"cov"]))
  }
  return(temp/s0)
}

tlist <- sort(unlist(lapply(simdatasam,function(x){
  tl <- length(x[,"time"])
  x[-tl,2]
} )))
est <- cumsum(sapply(tlist,function(t) lamesttest(t,1,simdatasam,100,10)))
plot(tlist,est,typ="l")
abline(0,1)
mean((est-tlist)^2)


##########################################################################################################
# Estimation for observation process
simdatares <- sapply(1:100,function(id) datasetgen(1,id))

kernelh <- function(u,h) ifelse(abs(u/h)<1,3/4*(1-(u/h)^2)/h,0)

S0_1_in <- function(t,gamma,data,h){
  if(t < data$censoring){
    weightsum <- sum(kernelh(t-data$obscov_times,h))
    if(weightsum==0){
      res <-0
    } else {
      res <- sum(kernelh(t-data$obscov_times,h)*as.vector(exp(gamma%*%data$covariates)))/
        sum(kernelh(t-data$obscov_times,h))
    }
  } else {
    res <-0
  }
  
  return(res)
}

S0_1 <- function(t,gamma,datacol,h){
  S0_1_v <- sapply(datacol, function(data) S0_1_in(t,gamma,data,h))
  return(sum(S0_1_v))
}


#S0_1(2,1,simdatasam,2)

S1_1_in <- function(t,gamma,data,h){
  if(t < data$censoring){
    weightsum <- sum(kernelh(t-data$obscov_times,h))
    if(weightsum==0){
      res <-0
    } else {
      res <- data$covariates %*% (kernelh(t-data$obscov_times,h)*as.vector(exp(gamma%*%data$covariates)))/
        sum(kernelh(t-data$obscov_times,h))
    }
  } else {
    res <-0
  }
  
  
  return(as.vector(res))
}

S1_1 <- function(t,gamma,datacol,h){
  S1_1_v <- matrix(sapply(datacol, function(data) S1_1_in(t,gamma,data,h)),nrow=length(gamma))
  return(rowSums(S1_1_v))
}

#S1_1(2,1,simdatasam,2)

ugamma_1 <- function(datacol,gamma,h){
  
}

Zsmooth <- function(data,t,h){
  weightsum <- sum(kernelh(t-data$obscov_times,h))
  if(weightsum==0){
    res <-0
  } else {
    res <- data$covariates %*% matrix(kernelh(t-data$obscov_times,h),nrow=1)/
      sum(kernelh(t-data$obscov_times,h))
  }
  return(res)
}


library("foreach")

ugamma_1 <-function(datacol,gamma,h){
  part1 <- foreach(data=simdatasam,.combine='+') %do% {
    foreach(time = data$meas_times, .combine='+') %do% {
      weightsum <- sum(kernelh(time-data$obscov_times,h))
      if(weightsum==0){
        res <-0
      } else {
        res <- data$covariates %*% kernelh(time-data$obscov_times,h)/weightsum
      }
      as.vector(res)
    }
  }
  
  tcol <- unlist(sapply(datacol, function(data) data$meas_times))
  part2 <- sum(sapply(tcol, function(t) S1_1(t,gamma,datacol,h)/S0_1(t,gamma,datacol,h)))
  return(part1-part2)
}


#================================================================================
S1_1 <- function(t,gamma,datacol,h){
  foreach(data=datacol,.combine="+") %do% {
    res <-0
    weightsum <- sum(kernelh(t-data$obscov_times,h))
    if(weightsum !=0){
      res <- data$covariates %*% as.vector(kernelh(t-data$obscov_times,h)*exp(gamma %*% data$covariates))/weightsum
    }
    as.vector(res)
  }
}

S0_1 <- function(t,gamma,datacol,h){
  foreach(data=datacol,.combine="+") %do% {
    res <-0
    weightsum <- sum(kernelh(t-data$obscov_times,h))
    if(weightsum !=0){
      res <- sum(kernelh(t-data$obscov_times,h)*exp(gamma %*% data$covariates))/weightsum
    }
    as.vector(res)
  }
}

ugamma_1 <-function(datacol,gamma,h){
  part1 <- foreach(data=datacol,.combine='+') %do% {
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
  
  tcol <- unlist(sapply(datacol, function(data) data$meas_times))
  part2 <- sum(sapply(tcol, function(t) S1_1(t,gamma,datacol,h)/S0_1(t,gamma,datacol,h)))
  return(part1-part2)
}

ugamma_1(simdatasam,1,2)

#=======================================================================================================

nbar <- sapply(meas_obs_list, function(meastime) sapply(meastime, function(t) mean(sapply(meas_obs_list, function(x) sum(t>x)))))



#plot(unlist(meas_obs_list),unlist(sapply(testlong$Xbar_list,function(x) x[,2])))

plot(loess.smooth(unlist(meas_obs_list),
                  unlist(response_list)- unlist(sapply(testlong$Xbar_list, function(x) x%*%  testlong[[1]])),
                  span=nsample^(-0.5),evaluation = 10000),typ="l",ylim=c(-2,2))
curve((sin(2 * pi * x)),add=T,col="red")

lines(loess.smooth(unlist(meas_obs_list),
                  unlist(response_list)- unlist(sapply(testlong$Xbar_list, function(x) x[,1] *  2+(1:length(x[,1])-1))),
                  span=nsample^(-0.5),evaluation = 10000),typ="l",col="green")
curve((sin(2 * pi * x)),add=T,col="red")

lines(loess.smooth(unlist(meas_obs_list),
                   unlist(response_list)- unlist(sapply(testlong$Xbar_list, function(x) x[,1] *  2))-
                     unlist(sapply(nbar, function(x) x*1)),
                   span=nsample^(-0.5),evaluation = 10000),typ="l",col="yellow")

lines(loess.smooth(unlist(meas_obs_list),
                  unlist(response_list)- unlist(sapply(testlong$Xbar_list, function(x) x%*%  c(2,1))),
                  span=nsample^(-0.5),evaluation = 10000),typ="l",col="blue")

lines(loess.smooth(unlist(meas_obs_list),
                   unlist(response_list)- unlist(sapply(testlong$Xbar_list, function(x) x[,1] *  testlong[[1]][1,1]))-
                     unlist(sapply(nbar, function(x) x*testlong[[1]][2,1])),
                   span=nsample^(-0.5),evaluation = 10000),typ="l",col="orange")

curve((sin(2 * pi * x)),add=T,col="red")

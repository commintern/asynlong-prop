

# Epanechnikov Kernel function

kernelh <- function(u, h)
  ifelse(abs(u / h) < 1, 3 / 4 * (1 - (u / h) ^ 2) / h, 0)

# Estimate for measurement process: gamma, lambda ----------------------------------------

gammaest <-
  function(kerMat,
           meas_obs_list,
           covar_list,
           censor_list,
           n,
           p) {
    # U gamma wrapper
    part1 <- ugamma1_C(kerMat, covar_list, n, p)
    ugamma <- function(gamma) {
      part1 -
        ugamma2_C(gamma, kerMat, meas_obs_list, covar_list, censor_list, n, p)
    }

    gammaestres <- nleqslv(rep(0, p), ugamma)
  }

# Estimation of Lambda


lambdaest <-
  function(gammahat,
           kerMat,
           meas_obs_list,
           covar_list,
           censor_list,
           n,
           p) {
    # U gamma wrapper
    dlambdaest <-
      dlambda_C(gammahat, kerMat, meas_obs_list, covar_list, censor_list, n, p)
    res <- cbind(meas_obs_list, lapply(dlambdaest, as.vector))
    return(res)
  }

# Estimation


asymptotic <-
  function(kerMat,
           meas_obs_list,
           covar_list,
           censor_list,
           response_list,
           thetahat,
           gammahat,
           n,
           p,
           h) {
    Xbar_res <-
      Xbar_c(
        gammahat,
        thetahat,
        kerMat = kerMat,
        meas_times = meas_obs_list,
        covariates = covar_list,
        censor = censor_list,
        n,
        p
      )

    #print(Xbar_res$XXtbar[[1]])

    #
    # S0_res <- Xbar_res$S0
    # Xbar_res <- Xbar_res$Xbar
    # Rs_res <- Rs_c(
    #   gammahat,
    #   thetahat,
    #   kerMat = kerMat,
    #   meas_times = meas_obs_list,
    #   covariates = covar_list,
    #   Xbar = Xbar_res,
    #   S0 = S0_res,
    #   censor = censor_list,
    #   n,
    #   p
    # )
    H_A_res <-
      H_A_c(
        gamma = gammahat,
        theta = thetahat,
        kerMat = kerMat,
        meas_times = meas_obs_list,
        covariates = covar_list,
        response = response_list,
        Xbar = Xbar_res$Xbar,
        XXtbar = Xbar_res$XXtbar,
        XXtZbar = Xbar_res$XXtZbar,
        S0 = Xbar_res$S0,
        censor = censor_list,
        n = n,
        p = p
      )

    #print(H_A_res)
    # Pmat <- P_D_res$P
    # Dmat <- P_D_res$D
    # DO Not forget the kernel part
    Asymest <-
      long_asy_c(
        gamma = gammahat,
        theta = thetahat,
        kerMat = kerMat,
        meas_times = meas_obs_list,
        covariates = covar_list,
        Xbar = Xbar_res$Xbar,
        XXtbar = Xbar_res$XXtbar,
        response = response_list,
        Hmat = H_A_res$Hmat,
        Amat = H_A_res$Amat,
        censor = censor_list,
        n = n,
        p = p
      )
    #Sigmaest <- Sigmaest
    # Bmat <-
    #   long_asy_B_c(
    #     gamma = gammahat,
    #     theta = thetahat,
    #     kerMat = kerMat,
    #     meas_times = meas_obs_list,
    #     covariates = covar_list,
    #     response = response_list,
    #     dlambda = dlambda_list,
    #     censor = censor_list,
    #     n = n,
    #     p = p
    #   )
    # Vtemp <- EEZV_c(
    #   gamma = gammahat,
    #   theta = thetahat,
    #   kerMat = kerMat,
    #   meas_times = meas_obs_list,
    #   covariates = covar_list,
    #   Xbar = Xbar_res$Xbar,
    #   XXtbar = Xbar_res$XXtbar,
    #   response = response_list,
    #   Hmat = H_A_res$Hmat,
    #   Amat = H_A_res$Amat,
    #   censor = censor_list,
    #   n = n,
    #   p = p
    # )
    vargammaest <- solve(H_A_res$Amat, Asymest$Vmat) %*% t(solve(H_A_res$Amat))
    #vargammaest <- solve(H_A_res$Amat, Vtemp) %*% t(solve(H_A_res$Amat))
    varest <- solve(Asymest$Bmat, Asymest$Sigma) %*% t(solve(Asymest$Bmat))
    return(rbind(as.vector(gammahat) + cbind(-1 * qnorm(0.975) * sqrt(diag(vargammaest)), qnorm(0.975) * sqrt(diag(vargammaest))),as.vector(thetahat) + cbind(-1 * qnorm(0.975) * sqrt(diag(varest)), qnorm(0.975) * sqrt(diag(varest)))))
  }

#kernel <- function(x,h) exp(-2*abs(x/h))/h
estasy <- function(dataset, kerFun, h, n, p) {
  # transform data
  covar_list <- lapply(dataset, function(x)
    x[["covariates"]])
  meas_obs_list <- lapply(dataset, function(x)
    x[["meas_times"]])
  obscov_times_list <-
    lapply(dataset, function(x)
      x[["obscov_times"]])
  censor_list <- sapply(dataset, function(x)
    x[["censoring"]])

  # R realization, slow
  #kerMat <- apply(expand.grid(1:n,1:n),1, function(ind) kernelh(outermin_C(dataset[[ind[2]]][[2]],dataset[[ind[1]]][[4]]),h))

  kerMat <- kerMatgen_C(meas_obs_list, obscov_times_list, h)

  #browser()
  gammaest_res <-
    gammaest(kerMat, meas_obs_list, covar_list, censor_list, n, p)
  # TODO consider add error handling and logging if not converges

  print(gammaest_res$x)

  lambdaest_res <-
    lambdaest(gammaest_res$x,
              kerMat,
              meas_obs_list,
              covar_list,
              censor_list,
              n,
              p)

  time_order <- order(unlist(lambdaest_res[, 1]))
  Lambdaest_res <-
    cbind(sort(unlist(lambdaest_res[, 1])), cumsum(unlist(lambdaest_res[, 2])[time_order]))

  lambda_persub <- numeric(length(time_order))
  lambda_persub[time_order] <- Lambdaest_res[,2]

  lambda_persub <- split(lambda_persub,rep(1:n,sapply(lambdaest_res[,1],length)))

  lambda_persub <- sapply(lambda_persub, function(x) diff(c(0,x)))

  response_list <- sapply(dataset, function(x)
    x[["Y"]])
  #kerMat <- kerMatgen_C(meas_obs_list,obscov_times_list,nsample ^ (-0.8))

  longest_res <- longest_c(
    gamma = gammaest_res$x,
    kerMat = kerMat,
    meas_times = meas_obs_list,
    covariates = covar_list,
    response = response_list,
    dlambda = NULL,
    censor = censor_list,
    n,
    p
  )

  CI_theta <- asymptotic(kerMat,
             meas_obs_list,
             covar_list,
             censor_list,
             response_list,
             longest_res[[1]],
             gammaest_res$x,
             lambda_persub,
             longest_res[[2]],
             n,
             p,h)

  #gmu0est <- cbind(unlist(meas_obs_list),unlist(response_list)-do.call(rbind,longest$Xbar_list) %*% longest[[1]])
  #gmu0est <- gmu0est[order(gmu0est[,1]),]
  rm(kerMat)
  return(
    list(
      gamma = gammaest_res$x,
      lambda0 = lambdaest_res,
      Lambda0 = Lambdaest_res,
      thetaest = longest_res[[1]],
      gmu0est = longest_res[[2]],
      CI_theta = CI_theta
    )
  )
}

#============ Purturbation ==========================

gammaest_pur <-
  function(kerMat,
           meas_obs_list,
           covar_list,
           censor_list,
           n,
           p,
           pur_weights) {
    # U gamma wrapper
    part1 <- ugamma1_pur_C(kerMat, covar_list, n, p, pur_weights)
    ugamma <- function(gamma) {
      part1 -
        ugamma2_pur_C(gamma,
                  kerMat,
                  meas_obs_list,
                  covar_list,
                  censor_list,
                  n,
                  p,
                  pur_weights)
    }

    gammaestres <- nleqslv(rep(0, p), ugamma)
  }


purtur_CI_one <-
  function(kerMat,
           meas_obs_list,
           covar_list,
           censor_list,
           response_list,
           n,
           p,
           pur_weights) {
    #browser()
    gammaest_res <-
      gammaest_pur(kerMat,
                   meas_obs_list,
                   covar_list,
                   censor_list,
                   n,
                   p,
                   pur_weights)

    longest_res <- longest_pur_c(
      gamma = gammaest_res$x,
      kerMat = kerMat,
      meas_times = meas_obs_list,
      covariates = covar_list,
      response = response_list,
      dlambda = NULL,
      censor = censor_list,
      n,
      p,
      pur_weights
    )

    return(c(gammaest_res$x, longest_res[[1]]))

  }

estasy_pur <- function(dataset, kerFun, h, n, p) {
  # transform data
  covar_list <- lapply(dataset, function(x)
    x[["covariates"]])
  meas_obs_list <- lapply(dataset, function(x)
    x[["meas_times"]])
  obscov_times_list <-
    lapply(dataset, function(x)
      x[["obscov_times"]])
  censor_list <- sapply(dataset, function(x)
    x[["censoring"]])

  response_list <- sapply(dataset, function(x)
    x[["Y"]])

  # R realization, slow
  #kerMat <- apply(expand.grid(1:n,1:n),1, function(ind) kernelh(outermin_C(dataset[[ind[2]]][[2]],dataset[[ind[1]]][[4]]),h))

  kerMat <- kerMatgen_C(meas_obs_list, obscov_times_list, h)

  point_est <- purtur_CI_one(kerMat,
             meas_obs_list,
             covar_list,
             censor_list,
             response_list,
             n,
             p,
             rep(1,n))

  CI_asym <- asymptotic(kerMat,
                        meas_obs_list,
                        covar_list,
                        censor_list,
                        response_list,
                        point_est[-(1:p)],
                        point_est[1:p],
                        n,
                        p,h)
  # CIrep <- apply(matrix(rexp(500*n,rate=1),ncol=n),1,function(ww) purtur_CI_one(kerMat,
  #                                                                meas_obs_list,
  #                                                                covar_list,
  #                                                                censor_list,
  #                                                                response_list,
  #                                                                n,
  #                                                                p,
  #                                                                ww))
  # varest <- var(t(CIrep))
  # CIres <- point_est + cbind(-1 * qnorm(0.975) * sqrt(diag(varest)),
  #                             qnorm(0.975) * sqrt(diag(varest)))

  CIres <- CI_asym



  rm(kerMat)
  return(
    list(
      point_est,
      CI_theta = CI_asym,
      CI_pur = CIres
    )
  )
}

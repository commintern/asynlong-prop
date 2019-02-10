broyden <- function(Ffun, x0, J0 = NULL, ...,
                    maxiter = 100, tol = .Machine$double.eps^(1/2)) {
  if (!is.numeric(x0))
    stop("Argument 'x0' must be a numeric (row or column) vector.")
  fun <- match.fun(Ffun)
  F <- function(x) fun(x, ...)
  y0 <- F(x0)
  if (length(x0) != length(y0))
    stop("Function 'F' must be 'square', i.e. from R^n to R^n .")

  # Compute once the Jacobian and its inverse
  # if (is.null(J0)) {
  #   #A0 <- jacobian(F, x0)
  # } else {
  #   A0 <- J0
  # }
  #
  # B0 <- inv(A0)
  # if (any(is.infinite(B0)))
  #   B0 <- diag(length(x0))
  B0 <- diag(length(x0))

  # Secant-like step in Broyden's method
  xnew <- x0 - B0 %*% y0 # cause B0 is I
  ynew <- F(xnew)

  s <- xnew - x0
  d <- ynew - y0
  B0 <- B0 + (d - B0 %*% s) %*% t(s) / c(t(s) %*% s)

  x0 <- xnew
  y0 <- ynew

  k <- 1
  while (k < maxiter) {
    #browser()
    print(as.vector(x0))
    print(sum(y0))

    if (norm(s, "F") < tol || norm(as.matrix(ynew), "F") < tol)  break

    # Update for next iteration
    #x0 <- xnew
    lambda <- 1
    p <- solve(B0,y0)
    xnew <- x0 - lambda * p
    ynew <- F(new)
    # Check descent to ensure globale convergence
    fnew <- sum(ynew^2)/2
    f0 <- sum(y0^2)/2
    if(fnew > f0){
      lambda2 <- lambda
      gp0 <- -f0*2
      lambda1 <- -gp0/(2*(fnew-g0-gp0))
      lambda1 <- max(0.1,lambda1)
      xnew <- x0 - lambda1 * p
      ynew <- F(new)
      fnew <- sum(ynew^2)/2
      f0 <- sum(y0^2)/2
      if(fnew > f0){

      }
    }

    #y0 <- ynew
    #ynew <- F(xnew)

    s <- xnew - x0
    d <- ynew - y0

    # Sherman-Morrison formula
    # B0 <- B0 + (s - B0 %*% d) %*% t(s) %*% B0 / c(t(s) %*% B0 %*% d)

    # Good Broyden update
    B0 <- B0 + (d - B0 %*% s) %*% t(s) / c(t(s) %*% s)



    k <- k + 1
  }

  if (k >= maxiter)
    warning(paste("Not converged: Max number of iterations reached."))

  fnew <- sqrt(sum(ynew^2))
  return(list(zero = c(xnew), fnorm = fnew, niter = k))
}

linesearch <- function(Ffun,F0,x0,B0){}

foo <- function(z) 1-z*exp(sum(z * c(1,2,3)))

#' @title Parameter estimation of Dirichlet-multinomial (DM) by Newton-Raphson (NR) algorithm.
#'
#' @description Given the data, we implemented NR to estimate parameters of DM which is faster than other methods.
#'
#' @param X a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the components.
#' @param init the initial value for the interative algorithm; if not given, we would like to use moment estimation.
#' @param iter the number of the interative steps.
#' @param conv the precison for stop the estimation.
#'
#' @return a list with the estimated parameters and log-likelihood value.

est.dm.NR <- function(X, init = NULL, iter = 1000, conv = 1e-6) {
  K <- ncol(X)
  if(any(colSums(X)==0)) {
    res <- list(alpha=rep(NA,K), loglik=-Inf)
    return(res)
  }
  X <- as.matrix(X[rowSums(X) >0, ])
  N <- nrow(X)
  rsx <- rowSums(X)

  if(is.null(init))  init <- mom.dm(X)
  init[is.na(init)] <- 1
  if(nrow(X) < 5)  init

  alpha0  <- init
  loglik0 <- loglik_dm(X, alpha0)

  result <- matrix(0, iter, K+1)
  result[1, ] <- c(alpha0, loglik0)
  for(i in 2:iter) {
    g <- N*digamma(sum(alpha0)) - sum(digamma(rsx + sum(alpha0))) +
      rowSums(digamma(t(X)+alpha0)) - N* digamma(alpha0)
    q <- rowSums(trigamma(t(X) + alpha0)) - N*trigamma(alpha0)
    z <- N* trigamma(sum(alpha0)) - sum(trigamma(rsx + sum(alpha0)))
    b <- sum(g/q)/(1/z + sum(1/q))
    hg <- (g-b)/q
    alpha1 <- alpha0 - hg

    if(any(alpha1<=0) | any(is.na(alpha1))) {result <- result[i-1, ]; break}

    loglik1 <- loglik_dm(X, alpha1)
    result[i, ] <- c(alpha1, loglik1)

    err <- (loglik1 - loglik0)/abs(loglik0+1e-10)
    if(err <= conv | is.na(err)) { result <- result[i, ]; break}

    alpha0  <- alpha1
    loglik0 <- loglik1
  }
  if(all(is.na(result))) result <- rep(1, K+1)
  # result
  res <- list(alpha=result[1:K], loglik=result[K+1])
  res
}

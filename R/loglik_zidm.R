#' @title log-likelihood function for Zero-inflated Dirichlet-multinomial (ZIDM) distribution.
#'
#' @description Given data X and the parameters, we could calcualte the log-likelihood value of them.
#'
#' @param X a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the components
#' @param alpha the corresponding parameters for Dirichlet-multinomial part.
#' @param pi the correpsonding parameters for Zero-inflated part
#'
#' @return the log-likelihood value for ZIDM


loglik_zidm <- function(X, alpha, pi) {
  # for binary, only (one pi) pi=0 is DM
  X <- as.matrix(X)
  X <- X[rowSums(X)>0,]
  N <- nrow(X)
  Z <- rep(0, N)
  Z[X[, 1]==0] <- 1
  rx <- rowSums(X)
  # lfdm <- N* lgamma(sum(alpha)) - sum(lgamma(rowSums(X)+sum(alpha))) + sum(lgamma(t(X)+alpha)) - N* sum(lgamma(alpha))
  lfdm1 <- lgamma(sum(alpha)) - lgamma(rx+sum(alpha)) + colSums(lgamma(t(X)+alpha)) -
    sum(lgamma(alpha)) + lgamma(rx+1) - rowSums(lgamma(X+1))
  res <- sum(log(pi*Z + (1-pi)*exp(lfdm1)))
  res
}

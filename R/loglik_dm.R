#' @title log-likelihood function for Dirichlet-multinomial (DM).
#'
#' @description Given data X and the parameters, we could calcualte the log-likelihood value of them.
#'
#' @param X a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the components.
#' @param alpha the corresponding parameters.
#'
#' @return the log-likelihood value for DM.

loglik_dm <- function(X, alpha) {
  N <- nrow(X)
  res <- N* lgamma(sum(alpha)) - sum(lgamma(rowSums(X)+sum(alpha))) +
    sum(lgamma(t(X)+alpha)) - N* sum(lgamma(alpha)) + sum(lgamma(rowSums(X)+1)) - sum(lgamma(X+1))
  return(res)
}

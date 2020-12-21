#' @title Zero-inflated Dirichlet (ZID) distribution.
#'
#' @description Generate ZID probability with given parameters, sample size, it will degenerate to zero-inflated beta with only two componenets.
#'
#' @param N a value for the sample size.
#' @param alpha a vector for the with K positive parameters.
#' @param pi a value or vector with length K-1, for the zero proportion of the first K-1 componenets.
#' @param seed a value for the seeds.
#' @return a probability matrix, with some zero.
#'
#' @examples
#' X <- rzidirichlet(N=10, alpha=rep(1, 5), pi=0.1)
#' X
#' @importFrom stats rbinom rbeta
#' @export
#'

rzidirichlet <- function(N=10, alpha, pi=0.1, seed=1) {
  # N, sample size
  # alpha, parameters for beta distribution
  # pi, zero proportion
  K <- length(alpha)-1
  if(length(pi) != K)  pi<- rep(pi[1], K) # zero propability equally on p-1 first variable and matching alpha
  if(any(alpha < 0)) stop('Alpha should be positive.')
  if(any(pi < 0) | any(pi>1)) stop('Pi should range between 0 and 1.')
  Pi <- matrix(0, N, K+1)
  for(i in 1:N) {
    z <- rep(0, K)
    set.seed(seed*i)
    delta <- rbinom(K, 1, pi)
    zi <- which(delta==0)
    if(length(zi)>0) z[zi] <- sapply(zi, function(x) { rbeta(1, alpha[x], sum(alpha[(x+1):(K+1)]), ncp = 0)})
    if(K==1) {
      Pi[i,] <- c(z, 1-z)
    }else{
      tpi <- sapply(2:K, function(x) z[x]*prod(1-z[1:(x-1)]))
      Pi[i,] <- c(z[1], tpi, 1-z[1]-sum(tpi))
    }
  }
  Pi
}

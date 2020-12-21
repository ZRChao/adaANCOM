#' @title Outlier detection
#'
#' @description To avoid abnormal value based on the transformed data of the log-ratios, we use the maximum of the log-ratios of non-zeros as the threshold.
#' @param X a data frame or matrix containing the count data, with two columns.
#' @param P a data frame or matrix containing the transformed data, with two columns.
#'
#' @return a list containing the index of the outliers and the value of the log-ratios both for noraml and abnormal value.
#' @examples
#'
#' # data generation
#' N <- 100
#' D1 <- D2 <- round(runif(N, 1, 100))
#' Pi <- rzidirichlet(N=N, alpha =  c(3, 6), pi=0.1)
#' which(Pi[,1]==0) # true zero
#' X1 <- matrix(NA, N, 2)
#' for(i in 1:N) X1[i, ] <- rmultinom(1, D1[i], prob=Pi[i, ])
#'
#' prop <- Smooth_X(X1, type='MIX', rel=TRUE, outlier=FALSE)
#'
#' res <- Outliers(X1, prop)
#' res
#' @export
#'


Outliers <- function(X, P) {

  X[X!=0] <- 1
  idx0 <- which(rowSums(X)==0)
  idx1 <- which(rowSums(X)==1)
  idx2 <- which(rowSums(X)==2)

  ratios <- log(P[,1]) - log(P[,2])

  r1 <- abs(ratios[idx1])
  r2 <- abs(ratios[idx2])

  liers <- c(idx0, idx1[r1>max(r2)])

  return(list(outliers=liers, ratios=ratios[-liers], abnormal=ratios[liers]))
}

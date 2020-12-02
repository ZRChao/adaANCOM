#' @title Log-ikelihood ratio test (lrt) between Multinomial (MN) and Dirichlet-mulitnomial (DM) distribution.
#'
#' @description Since MN nested in DM, we construct lrt between MN and DM, the degree of freedom equals to 1.
#'
#' @param X a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the components.
#' @param fdm the log-likelihood value of DM.
#'
#' @return p-value for the test.
#' @import stats

LiK_MN <- function(X,  fdm=NULL) { # mn vs dm
  ry <- rowSums(X)
  X  <- as.matrix(X[ry>1, ])
  ry <- ry[ry>1]
  pk <- colSums(X)/sum(ry)
  fmn <- sum(lgamma(ry+1)) - sum(lgamma(X+1)) + sum(t(X)*log(pk))
  if(is.null(fdm)) fdm <- est.dm.NR(X)[[2]]
  p2 <- pchisq(-2*(fmn - fdm),   df=1, lower.tail = F)
  if(is.na(p2)) p2 <- 1
  p2
}

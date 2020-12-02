#' @title Moment estimation for Dirichlet-multinomial
#'
#' @description Moment estimation for Dirichlet-multinomial.
#'
#' @param X a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the components
#'
#' @return a estimated vector
#' @export
#'


mom.dm <- function(X) {
  # row is sample
  # delete all zero sample
  X <- as.matrix(X); rX <- rowSums(X)
  X <- X[rX > 0, ];  tx<- X/rowSums(X)
  p1<- tx[,1]
  p1 <- p1[!is.na(p1)]
  init <- (mean(p1) - mean(p1^2))/(var(p1))*colMeans(tx, na.rm = T)
  init[init<=0] <- 1
  init
}

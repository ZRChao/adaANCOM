#' @title Log-ikelihood ratio test (lrt) between Dirichlet-mulitnomial (DM) and Zero-inflate DM (ZIDM) distribution.
#'
#' @description Since DM nested in ZIDM, we construct lrt between them, the degree of freedom equals to the number of the conmponent minus 1.
#'
#' @param X a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the components.
#' @param fdm the loglikelihood value of DM.
#' @param  fzidm the loglikelihood value of ZIDM.
#'
#' @return p-value for the test.
#'
#' @export
#'

LiK_DM <- function(X, fdm=NULL, fzidm=NULL) { # dm vs zidm
  ry <- rowSums(X)
  X  <- as.matrix(X[ry>1, ])
  p <- ncol(X)
  if(is.null(fzidm)) fzidm <- est.zidm.EM(X)[[3]]
  if(is.null(fdm)) fdm <- est.dm.NR(X)[[2]]
  p1 <- pchisq(-2*(fdm - fzidm), df=p-1, lower.tail = F)
  if(is.na(p1)) p1 <- 1
  p1
}

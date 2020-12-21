#' @title Log-ikelihood ratio test (lrt) between Multinomial, Dirichlet-mulitnomial (DM), and Zero-inflate DM (ZIDM) distribution.
#'
#' @description Since MN < DM < ZIDM,
#'
#' @param X a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the components.
#' @note also see \link{lrt_DM2ZIDM} or \link{lrt_MN2DM}.
#' @return p-value for the test and the log-likelihood value
#' @export
#' @import stats
#'

LikRTest <- function(X) {
  Y <- X
  p <- ncol(X)
  ry <- rowSums(Y)
  fzidm=NULL; fdm=NULL
  Y  <- as.matrix(Y[ry>1, ])
  if(is.null(fzidm)) {
    ea <- est.zidm.EM(Y); #a <- ea[[1]]; pi <- ea[[2]]
    fzidm <- ea[[3]]
  }
  pk <- colSums(Y)/sum(ry)
  fmn <- sum(lgamma(ry+1)) - sum(lgamma(Y+1)) + sum(t(Y)*log(pk))
  #fmn <- sum(lgamma(ry+1)) - sum(lgamma(Y+1)) + sum(Y*log(Y/ry), na.rm = T)
  if(is.null(fdm)) fdm <- est.dm.NR(Y)[[2]]
  # fzidm <- loglik_zidm(Y, a, pi=pi)

  p1 <- pchisq(-2*(fdm - fzidm), df=p-1, lower.tail = F)
  p2 <- ifelse(is.na(fmn), 1, pchisq(-2*(fmn - fdm),   df=1, lower.tail = F))
  # p3 <- pchisq(-2*(fmn - fzidm), df=2, lower.tail = F)
  # res <- data.frame(Test=c('H0:DM-H1:ZIDM', 'H0:MN-H1:DM'), pvalue=c(p1, p2))
  res <- list(pvalue=c(p1, p2), loglik=c(fmn, fdm, fzidm))
  # res <- c(p1, p2)
  return(res)
}

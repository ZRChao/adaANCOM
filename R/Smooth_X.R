#' @title Zero smooth.
#'
#' @description We used the posterior mean transfromation to smooth zero. If the prior is not given,
#' we could use likelihood ratio test to chooce the best fitted model and the corresponding posterior mean to smooth zeros.
#'
#' @param X a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the components
#' @param type the smoothing methods which could be specified including 0.5, DM, ZIDM, or MIX which using lrt to choose the best fitted one.
#' @param rel If true, the output will be the posterior; otherwise, we product the posterior by the sample depth to output the smoothed count.
#' @param outlier If true, remove the outliers, otherwise keep them.
#'
#' @return smoothed data matrix
#' @export
#'


Smooth_X <- function(X, type='DM', rel=F, outlier=F) {
  zi<- which(rowSums(X)>0)
  P <- matrix(NA, nrow(X), ncol(X))
  X <- as.matrix(X[zi, ])
  if(length(zi)<5){ # for little sample size estimation with big variance
    Y <- X + 0.5
    P[zi, ] <- Y
    if(rel)  {
      if(length(zi) == 1) sy <- sum(Y)
      if(length(zi) > 1) sy <- rowSums(Y)
      P[zi, ] <- Y/sy
    }
    return(P)
  }
  a <- NULL

  zero_detect <- function(X, a, pi, sig=0.5) {
    # idx for non zero
    ix <- which(X[,1]==0)
    tmp <- exp(lbeta(a[1], a[2]+X[ix,2]) - lbeta(a[1],a[2]))
    px <- pi/(pi+(1-pi)*tmp)
    if(nrow(X)-length(ix)<5) {
      sig <- 0.5
    }else{
      px1 <- (X[-ix,1]+a[1])/(rowSums(X[-ix,])+sum(a))
      sig <- sort(px1, decreasing = T)[2]
    }
    ix[px>sig]
  }

  if(type == 'MIX')  {
    em <- est.zidm.EM(X)
    nr <- est.dm.NR(X)
    p1 <- LiK_DM(X, fdm=nr[[2]], fzidm=em[[3]])
    p2 <- LiK_MN(X, fdm=nr[[2]])
    pi <- em[[2]]
    if(p1<0.05){
      type <- 'ZIDM'; a <- em[[1]]
      if(p2>0.05 | sum(a)>99) type <- a <- 0.5
    }else{
      type <- 'DM'; a <- nr[[1]]
      if(p2>0.05) type <- a <- 0.5
    }
  }

  if(type == 'DM')  {
    if(is.null(a)) nr <- est.dm.NR(X)
    a <- nr[[1]]
    p2 <- LiK_MN(X, fdm=nr[[2]])
    if(sum(a)>99 | p2>0.05)  a <- 0.5
    if(sum(a) < 1) a <- 0.5
    Y <- t(t(X)) + a
    Y[rowSums(X)< sum(a)+1, ] <-  X[rowSums(X)< sum(a)+1, ] + 0.5 #X[rowSums(X)< sum(a)+1, ]
  }
  if(is.numeric(type))  {
    Y <- t(t(X)) + type
  }

  if(type == 'ZIDM') {
    Y <- X
    X1 <- X[X[, 1]==0, 2]
    t1 <- lbeta(a[1]+1, a[2]+X1) - lbeta(a[1],a[2])
    t2 <- lbeta(a[1], a[2]+X1) - lbeta(a[1],a[2])
    Y[X[,1]==0, 1] <- (1-pi)*exp(t1)/(pi+(1-pi)*exp(t2))*(X1+sum(a))

    Y[X[,1]!=0, 1] <- X[X[, 1]!=0, 1]+a[1] # /(rowSums(X2) + sum(a))
    Y[, 2] <- rowSums(X) + sum(a) - Y[, 1]

    # idx <- which(X[, 1]==0 & rowSums(X)<sum(a)+1)
    Y[rowSums(X) < sum(a)+1, ] <- X[rowSums(X)<sum(a)+1, ]
    if(outlier) {
      zo <- zero_detect(X, a=a, pi=pi)
      Y[zo, ] <- X[zo, ]
    }
    # P[zi, 1] <- c(Y1, Y2); P[zi, 2] <- 1 - P[zi, 1]
  }
  P[zi, ] <- Y
  if(rel) P[zi, ] <- Y/rowSums(Y)
  P
}

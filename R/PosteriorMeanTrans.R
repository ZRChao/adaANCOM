#' @title Posterior mean transformation
#'
#' @description Data transformation based on the posterior mean.
#' @param otu the OTU count table.
#' @param tree the phylogenetic tree.
#' @param model smooth methods, also see \link{Smooth_X}.
#'
#' @return Returns two lists:
#' \itemize{
#' \item count count abundance for internal nodes V, leaf nodes L and count on each subtree Y.
#' \item probability abundance for internal nodes V, leaf nodes L and count on each subtree Y.
#' }
#'
#' @importFrom ape is.binary
#' @export
#'
#'


PosteriorMeanTrans <- function(otu, tree, model='MIX') {
  if(!is.binary(tree)) stop('Binary tree only !')
  p <- ncol(otu)
  if(is.null(tree$Nnode)) tree$Nnode <- p-1
  Nnod <- tree$Nnode
  taxa <- colnames(otu)
  tip <- tree$tip.label
  edge <- tree$edge
  if(p != length(tip)) stop('OTU table not matched to the tree!')
  if(is.null(taxa)) taxa <- tip
  if(length(intersect(taxa, tip))!=p) stop('OTU table not matched to the tree!')

  Smooth_XX <- function(X, type='DM', sig=0.01) {
    X <- as.matrix(X)
    cx<- colSums(X)
    P <- matrix(0, nrow(X), ncol(X))
    zi<- which(rowSums(X)>0)

    if(cx[1]==0 & cx[2]==0) {return(P)}
    if(cx[1]==0 | cx[2]==0 | length(zi)<5) { # for little sample size estimation with big variance
      Y <- X + 0.5
      P <- Y/rowSums(Y)
      P[-zi, ] <- 0
      return(P)
    }

    X <- as.matrix(X[zi, ])
    a <- NULL
    Ts <- F
    alpha <- sig

    if(type == 'MIX')  {
      em <- est.zidm.EM(X)
      em2 <- est.zidm.EM(X[, 2:1])
      nr <- est.dm.NR(X)
      if(is.na(em$loglik)) em <- nr
      if(is.na(em2$loglik)) em2 <- nr
      if(em2$loglik > em$loglik)  {em <- em2; Ts <- T}
      if(Ts) X <- X[, 2:1]

      p1 <- lrt_DM2ZIDM(X, fdm=nr[[2]], fzidm=em[[3]])
      p2 <- lrt_MN2DM(X, fdm=nr[[2]])
      pi <- em[[2]]
      if(p1<alpha){
        type <- 'ZIBB'; a <- em[[1]]
        if(p2>alpha | sum(a)>99) type <- a <- 0.5
      }else{
        type <- 'BB'; a <- nr[[1]]
        if(p2>alpha) type <- a <- 0.5
      }
    }

    if(type == 'BB')  {
      if(is.null(a)) nr <- est.dm.NR(X)
      a <- nr[[1]]
      p2 <- lrt_MN2DM(X, fdm=nr[[2]])
      if(sum(a)>99 | p2>0.05)  a <- 0.5
      if(sum(a) < 1) a <- 0.5
      Y <- t(t(X)) + a
      Y[rowSums(X)< sum(a)+1, ] <-  X[rowSums(X)< sum(a)+1, ] + 0.5 #X[rowSums(X)< sum(a)+1, ]
    }
    if(is.numeric(type))  {
      Y <- t(t(X)) + type
    }

    if(type == 'ZIBB') {
      Y <- X
      idx <- which(X[,1]==0)
      X1 <- X[idx, 2]
      t1 <- lbeta(a[1]+1, a[2]+X1) - lbeta(a[1],a[2])
      t2 <- lbeta(a[1], a[2]+X1) - lbeta(a[1],a[2])
      Y[X[,1]==0, 1] <- (1-pi)*exp(t1)/(pi+(1-pi)*exp(t2))*(X1+sum(a))

      Y[-idx, 1] <- X[-idx, 1]+a[1] # /(rowSums(X2) + sum(a))
      Y[, 2] <- rowSums(X) + sum(a) - Y[, 1]

      Y[rowSums(X) < sum(a)+1, ] <- X[rowSums(X)<sum(a)+1, ]+0.5
    }
    if(Ts) Y <- Y[, 2:1]
    P[zi, ] <- Y/rowSums(Y)
    return(P)
  }

  AccumProb <- function(tree, b){
    edge <- tree$edge
    n <- nrow(edge)
    if(length(b) != n) stop('Must same length as edges !\n')
    node <- unique(edge[,1])
    p <- length(tree$tip.label)
    for(i in 1:length(node)) {
      bi <- b[edge[, 1]==(i+p)]
      if(any(bi<0)) {stop('The probability should not be negative !')}
      if(sum(bi)!=0) {
        if(abs(sum(bi) -1) > 1e-6) {stop('The probability should sum to 1.')}
      }
    }

    pm <- bm <- b
    for(i in 1:n){
      idx <- which(edge[i, 1] == edge[, 2])
      if(length(idx)==0) {
        pm[i] <- bm[i]
      }else{
        pm[i] <- ifelse(bm[i]==0, pm[idx]/2, bm[i]*pm[idx])
      }
    }
    if(abs(sum(pm[edge[,2]<= p]) - 1)>1e-8) warning('Accumulative probability summation not equals to 1 for the leaves.')
    pm
  }

  N <- nrow(otu)
  # count
  Yotu <- Ytreefun(otu, tree)
  Y_count <- Y_prop <- Yotu[[1]]
  V_count <- Yotu[[2]][, (p+1):(p+Nnod)]

  # probability
  Altb <- matrix(NA, N, p+Nnod)
  Altb[, p+1] <- 1
  for(i in 1:Nnod) {
    Y_prop[[i]] <- tmp <- Smooth_XX(X = as.matrix(Y_count[[i]]), type=model)
    Altb[, edge[edge[,1]==(i+p),2]] <- tmp
  }

  Alltab0 <- Altb[, edge[,2]]
  for(i in 1:N) Alltab0[i, ] <- AccumProb(tree, Alltab0[i, ])
  V_prop <- Alltab0[, edge[,2]>p]
  L_prop <- Alltab0[, edge[,2]<(p+1)]

  res <- list(count=list(L=otu, V=V_count, Y=Y_count),
              probability=list(L=L_prop, V=V_prop, Y=Y_prop))
  return(res)
}

#' @title Zero smooth on the tree.
#'
#' @description Smooth on each internal node indenpendently.
#' @param Yotu a list which containing all data matrix on each internal nodes.
#' @param tree the phylogenetic tree.
#' @param type smooth methods, also see \link{Smooth_X}.
#' @param rel If true, the output will be the posterior; otherwise, we product the posterior by the sample depth to output the smoothed count.
#' @param outlier If true, remove the outliers, otherwise keep them.
#'
#' @return a list, one is the smoothed data on each internal nodes and other one is a matrix which consist of all internal and leaf nodes.
#'
#'

Smooth_Tree <- function(Yotu, tree, type='DM', rel=F, outlier=F) {
  Nnod <- tree$Nnode
  p <- length(tree$tip.label)
  edge <- tree$edge
  if(length(Yotu) != Nnod) stop('length different between list and taxa !')
  Altb <- matrix(NA, nrow(Yotu[[1]]), p+Nnod)
  colnames(Altb) <- as.character(1:(p+Nnod))
  for(i in 1:Nnod) {
    Y <- as.matrix(Yotu[[i]])
    X <- Smooth_X(X=Y, type=type, rel=rel, outlier = outlier)
    childs <- as.character(edge[edge[, 1] == (i+p),2])
    Altb[, childs] <- Yotu[[i]] <- X
  }
  res <- list(Ytree=Yotu, Altb=Altb)
  res
}

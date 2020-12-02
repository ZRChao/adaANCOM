#' @title Transform the count on each internal node to the leafs.
#'
#' @description Given a list which consist of all internal counts, we extract the count of the leaves.
#'
#' @param Yotu a list, containing each count matrix on the subtree.
#' @param tree the phylogenetic tree.
#' @return a OTU count table on the leafs.
#' @note also see \link{Ytreefun}.
#'
#' @export
#'

Ytree2otu <- function(Yotu, tree) {
  Nnod <- tree$Nnode
  p <- length(tree$tip.label)
  edge <- tree$edge
  if(length(Yotu) != Nnod) stop('length different between list and taxa !')
  Altb <- matrix(NA, nrow(Yotu[[1]]$count), p+Nnod)
  colnames(Altb) <- as.character(1:(p+Nnod))
  for(i in 1:Nnod) {
    Y <- as.matrix(Yotu[[i]]$count)
    childs <- as.character(edge[edge[, 1] == (i+p),2])
    Altb[, childs] <- Yotu[[i]] <- Y
  }
  otu <-Altb[, 1:p]
  otu
}

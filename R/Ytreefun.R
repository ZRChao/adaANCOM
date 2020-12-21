#' @title Calculate the count on each subtree.
#'
#' @description Follow the additive foumla, each internal counts equals to summation of his children's.
#'
#' @param otu.tab a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.
#' @param tree the phylogenetic tree.
#' @return Returns a list, containing each count matrix on the subtree.
#' @examples
#' library(phyloseq)
#' data(COMBO)
#' otu.tab <- t(otu_table(COMBO))
#' tree <- phy_tree(COMBO)
#' ytree <- Ytreefun(otu.tab, tree)
#' length(ytree)
#' @export
#' @import phyloseq
#'
#'


Ytreefun <- function(otu.tab, tree) {
  # calculate counts on each subtree
  edge <- tree$edge
  otus <- as.character(tree$tip.label);
  n.otu <- ncol(otu.tab)
  N <- nrow(otu.tab)
  if(is.null(otus)) otus <- tree$tip.label <- as.character(1:n.otu)
  otu1 <- colnames(otu.tab)
  if(is.null(otu1)) otu1 <- colnames(otu.tab) <- as.character(1:n.otu)
  if(is.null(tree$node.label)) taxa <- tree$node.label <- as.character((n.otu+1):(n.otu+tree$Nnode))
  if(length(intersect(otu1, otus)) != n.otu) stop('Table and Tree tips not mathced ! ')
  otu.tab <- as.matrix(otu.tab[, otus])

  node_order <- sapply(1:tree$Nnode, function(x) sum(find_descent(tree$edge, x+n.otu)<=n.otu))
  node_order <- order(node_order) + n.otu
  all_table  <- matrix(NA, N, n.otu+tree$Nnode)
  all_table[, 1:c(n.otu+1)] <- cbind(otu.tab, rowSums(otu.tab))
  ytree <- list()
  for(n in node_order) {
    chi <- tree$edge[tree$edge[,1]==n,2]
    ytree[[n-n.otu]] <- y <- all_table[, chi]
    all_table[, n] <- ry <- rowSums(y)
  }
  names(ytree) <- tree$node.label
  res <- list(Ytree=ytree, allcount=all_table)
  res
}

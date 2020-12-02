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
#' @import phyloseq
#' @export


Ytreefun <- function(otu.tab, tree) {
  # calculate counts on each subtree
  find_descent <- function(edge, v) {
    # find all the descents
    childs <- c(); par <- v
    while(length(par) >0){
      ch <- unlist(sapply(par, function(x) edge[edge[,1]==x,2]))
      par <- ch[ch>min(edge[,1])]
      childs <- c(childs, ch)
    }
    childs
  }

  n.otu  <- length(tree$tip.label)
  node_order <- sapply(1:tree$Nnode, function(x) sum(find_descent(tree$edge, x+n.otu)<=n.otu))
  node_order <- order(node_order) + n.otu
  all_table  <- matrix(NA, nrow(otu.tab), n.otu+tree$Nnode)
  all_table[, 1:c(n.otu+1)] <- cbind(as.matrix(otu.tab[, tree$tip.label]), rowSums(otu.tab))
  ytree <- list()
  for(n in node_order) {
    chi <- tree$edge[tree$edge[,1]==n,2]
    ytree[[n-n.otu]] <- y <- all_table[, chi]
    all_table[, n] <- ry <- rowSums(y)
  }
  names(ytree) <- tree$node.label
  ytree
}

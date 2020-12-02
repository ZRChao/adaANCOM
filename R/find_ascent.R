#' @title  Find ascentors
#'
#' @description Given the edge from phylo object and the number of the node, return its all ascentors
#' @param edge the edge from phylo object
#' @param v the number of the node
#' @return all ascentors
#' @examples
#' library(phyloseq)
#' data(COMBO)
#' tree <- phy_tree(COMBO)
#' find_ascent(tree$edge, 1)
#' @export
#'



find_ascent <- function(edge, v) {
  # find the ascents including v at first
  ascent <- v
  root <- edge[1,1]
  while(root != v) {
    par <- edge[edge[, 2] == v, 1]
    ascent <- c(ascent, par)
    v <- par
  }
  ascent
}

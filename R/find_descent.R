#' @title  Find descent
#'
#' @description Given the edge from phylo object and the number of the internal node, return its all descents.
#' @param edge the edge from phylo object
#' @param v the number of the internal node
#' @return all children nodes
#' @examples
#' library(phyloseq)
#' data(COMBO)
#' tree <- phy_tree(COMBO)
#' find_descent(tree$edge, 55)
#' @export
#' @import phyloseq

find_descent <- function(edge, v) {
  # find all the descents
  childs <- c()
  root <- min(edge[,1])
  if(v<root){
    childs <- v # return himself
  }else{
    par <- v
    while(length(par) >0){
      ch <- unlist(sapply(par, function(x) edge[edge[,1]==x,2]))
      par <- ch[ch>root]
      childs <- c(childs, ch)
    }
  }
  childs
}

#' @title Example dataset
#'
#' @description A example data with a binary tree with 10 leaves,
#' we assumed each internal nodes follow binomial distribution and
#' the probabilities for each edge store in the parameters with node 14
#' to his children be differential between two groups, one is (0.5, 0.5) while other groups is (0.55, 0.45).
#'Then the first three leaves are non-differential while left are differential.
#'
#' @docType data
#'
#' @usage data(example_data)
#'
#' @format A list contain an object of class \code{"phyloseq"} and a data.frame for the corresponding parameters.
#'
#' @keywords datasets
#'
#'
#'
#' @examples
#' data(example)
#' example_data$data
#'


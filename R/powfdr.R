#' @title Evaluation of the model.
#'
#' @description Given the statistics and threshold, we get the positive results.
#' Then compared to the true signal, we could calculate the recall, precision, and F1 score.
#'
#' @param p.adj the adjusted p-value or statistics.
#' @param dif the true differential conmonent.
#' @param a the significant level.
#'
#' @return a vector consist of recall, precision, and F1 score squentially.
#' @export
#'

powfdr = function(p.adj, dif = NULL, a = 0.05){
  reject = which(p.adj <= a)
  TP = length(intersect(reject, dif))
  FP = length(reject) - TP # type I error
  FN = length(dif) - TP    # type II error
  TN = length(p.adj) - TP-FP-FN

  recall = TP/(TP+FN) #power,sensitivity, true positive rate
  fdr = ifelse(length(reject) == 0, 0, FP/(TP+FP))
  fpr <- FP/(FP+TN)  # false positive rate
  spc <- 1-fpr  # specificity
  precision = ifelse(length(reject) == 0, 1, TP/(TP+FP)) # 1-fdr
  acc = (TP+TN)/length(p.adj)
  f1 <- ifelse((recall+precision)==0, 0, 2*(recall*precision)/(precision+recall))
  return(c(recall, precision, f1))
}

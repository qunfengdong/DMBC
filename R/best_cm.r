#' @title  Confusion Matrix for the Best Model
#'
#' @description Given an \code{loocv()} output and \code{Cal_AUC()} output, this function calculates the confusion matrix for the label produced from the model with optimized number of features.
#' @param CV Output object from loocv().
#' @param auc_out Output object form Cal_AUC()
#' @export
#' @return A list of confusionMatrix class.
#'
#' @examples
#' data(training)
#' cv = loocv(training)
#' auc_out = Cal_AUC(cv)
#' best_cm(CV=cv,auc_out=auc_out)
#'
best_cm <- function(CV=cv, auc_out=auc_out){
  bestrank <- auc_out$NumberOfFeature[which.max(auc_out$AUC_Area)]
  tmp <-lapply(cv,function(x) x[x$feature_rank==bestrank,])
  file_rank = data.frame(Reduce(function(...) merge(...,all=T),tmp))

  pred <- as.numeric(as.numeric(as.vector(file_rank$Type1_Posterior_Prb)) > as.numeric(as.vector(file_rank$Type2_Postereior_Prb)))
  truth <- as.numeric(as.vector(file_rank$Type1Label))

  confusionMatrix(pred,truth)

}

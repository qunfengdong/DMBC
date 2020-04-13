#' Example Training Data
#'
#' This is a simulated testing dataset with label for each sample.Relative abundance is used in DMBC version 1.1.0
#'data(training)
#'data(test)
#'CV=loocv(training)
#'auc_out =Cal_AUC(CV)
#'dmbc_predict(data=Meta_data,testset=test,auc_out=auc_out) ##delete
#'dmbc_predict(data=training,testSet=test,auc_out=auc_out)
#'
#' @docType data
#' @usage data(training)
#'
#' @keywords datasets
#' @examples
#' data(training)

"training"

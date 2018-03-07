#' @title 10-fold cross-validation
#' @description Implementation of 10-fold cross-validataion. It takes in an entire dataset, with first column the sample IDs, second column the group variable (classificatoin variable/disease type), using "Type" as its colum names, and third column the beginning of count table.
#' @param data A data frame validation set with first column the sample IDs, second column the group variable (classificatoin variable/disease type), and third column the beginning of count table.
#' @param type_col An index indicating at which column is group/type variable. The default is the 3rd column.
#' @param col_start An index indicating at which column is the beginning of bacteria (features) data. Default is the 2nd column.
#' @param Cutoff_mean The minimum average relative abundance allowed in filtering step. Default is 0.0005.
#' @param Cutoff_ratio The non-zero ratio cutoff in filtering features.Default value is 0.1.
#' @param totalReadsCutoff The minimum allowed total reads per sample. Any sample has less than this number of total reads will be removed. Default is 500.
#' @param Cutoff_pvalue The maximal P value allowed for a given feature to be remained in the list of selected features.
#' @export
#' @return A list of CalPrb() results. This output will be further used for Cal_AUC() in order to estimate the final model.
#' @examples
#' data(training)
#' tfcv(training)

tfcv <- function(data=data, type_col=2, col_start=3,Cutoff_mean=0.0005,Cutoff_ratio=0.1,totalReadsCutoff=500, Cutoff_pvalue = 0.5){
  res <- list()
  flds <- createFolds(data[,type_col],10,returnTrain = TRUE)
  for (i in 1:length(flds)){
    testSet <- data[-flds[[i]],]
    training <- data[flds[[i]],]
    FS_out <- FS(training=training,type_col=type_col,col_start=col_start,Cutoff_mean=Cutoff_mean,Cutoff_ratio=Cutoff_ratio,totalReadsCutoff=totalReadsCutoff, Cutoff_pvalue = Cutoff_pvalue)
    res[[i]] <- CalPrb(FS_out=FS_out,testSet=testSet,col_start=col_start,type_col=type_col,HighestRank=nrow(FS_out$Feature))
  }
  return(res)
}

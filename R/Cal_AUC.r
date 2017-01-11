#' @title Calculate Area Under the ROC Curve from Leave-One-Out Cross-Validation
#'
#' @description
#' This function calculates the AUC for each selected subset of features from DMBC cross-validation results, and output the set of features which provide the greatest AUC.
#' @param CV The output object from loocv()
#' @export
#' @return A data frame of 6 columns, comparison group1 label, comparison group2 label, number of features included in the model, AUC, AUC * prior probability, and selected features.
#' @examples
#' data(training)
#' Cal_AUC(loocv(training))


Cal_AUC <- function(CV=cv){
  HighestRank_allRow <- min(sapply(CV,nrow))
  auc_res <- list()
  for (NumberOfFeature in 1:HighestRank_allRow) {
    tryCatch({
      tmp <- lapply(CV,function(x) data.frame(x[x$feature_rank == NumberOfFeature,]))
      RankData = Reduce(function(...) merge(...,all=T),tmp)
     # AUC_Area= auc(roc(as.numeric(as.vector(RankData$Poster_Type1)), as.factor(RankData$Type1Label)))
      #AUC_Prior_Area= auc(roc(as.numeric(as.vector(RankData$Poster_Prio_Type1)), as.factor(RankData$Type1Label)))
      AUC_Area= auc(roc(as.numeric(as.vector(RankData$Type1_Posterior_Prb)), as.factor(RankData$Type1Label)))
      fnum <- sort(table(unlist(strsplit(as.character(RankData$SelectedFeatures),";"))),decreasing = T)
      features <- names(fnum[1:NumberOfFeature])
     # auc_res[[NumberOfFeature]] <- data.frame(RankData[1,1:2], NumberOfFeature, AUC_Area, AUC_Prior_Area,Features=paste(features,collapse = ";"))
      auc_res[[NumberOfFeature]] <- data.frame(RankData[1,1:2], NumberOfFeature, AUC_Area,Features=paste(features,collapse = ";"))

    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  Reduce(function(...) merge(...,all=T),auc_res)
}

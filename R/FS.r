#' @title Filter and Feature Selection based on Wilcoxon Rank-Sum test
#'
#' @description
#' Given a training set, this function performs feature selection based on several thresholds: (1). Average relative abudance in each cohort class (minimum relative abundance by default is 0.25\%), (2). Total reads per sample (minimum reads per sample is 500 by default), (3). Non-zero ratio out of all samples (By default, at least 10\% of the samples should have have non zero value.)
#' @param training A data frame of training set.
#' @param type_col An index indicating at which column is group/type variable. The default is the 3rd column.
#' @param col_start An index indicating at which column is the beginning of bacteria (features) data. Default is the 2nd column.
#' @param Cutoff_mean The minimum average relative abundance allowed in filtering step. Default is 0.0005.
#' @param Cutoff_ratio The non-zero ratio cutoff in filtering features.Default value is 0.1.
#' @param totalReadsCutoff The minimum allowed total reads per sample. Any sample has less than this number of total reads will be removed. Default is 500.
#' @param Cutoff_pvalue The maximum Pvalue allowed for a givien feature to be remained on the list of the selected features.
#' @return A list of 2: Feature and CountData.
#' @return \code{Feature}       A list of selected features sorted by their Wilcoxon P values.
#' @return \code{CountData}     A data frame containing balanced data.
#' @examples
#' data(training)
#' FS(training)
#' @export

#### use relative abundance *10000 do feature selection and following test

FS <- function(training=data,type_col=2,col_start=3,Cutoff_mean=0.001,Cutoff_ratio=0.1,totalReadsCutoff=500, Cutoff_pvalue = 1){
  ### Filter out samples less than totalReadsCutoff
  col_end=ncol(training)

  Row_freq = round(training[,col_start:col_end])
  RowFreqMeta = cbind(training[,1:(col_start-1)],Row_freq)
  fold = 10000/100

  taxalist <- list()
  for(taxa in col_start:col_end) {
    RowFreqMeta[,type_col] = droplevels(RowFreqMeta[,type_col])
    GRP <- split(RowFreqMeta[,taxa],RowFreqMeta[,type_col])
    GRP_abs <- split(training[,taxa],RowFreqMeta[,type_col])

    Mean = sapply(GRP,mean)
    MeanType1 = Mean[1]/fold
    MeanType2 = Mean[2]/fold
    #extract mean and SD for the taxa in the group 1 and group 2
    SD = sapply(GRP,sd)
    SDType1 = SD[1]/fold
    SDType2 = SD[2]/fold

    #extract all the taxa corresponding to group1 or group2
    group1 = GRP_abs[[1]]
    group2 = GRP_abs[[2]]

    Disease = levels(droplevels(training[,type_col]))

    non_zero1 = sum( group1!= 0)
    non_zero2 = sum( group2!= 0)
    ratio1 = non_zero1/length(group1)
    ratio2 = non_zero2/length(group2)

    #if all the sample has 0 in one taxa, fake the count of one sample into 1, for the Dirichlet purpose later. Otherwise, the Dirichlet-Multinomial distribution estimation gets weird
    if(length(group1) == length(group1[group1<1])) {
      RowFreqMeta[RowFreqMeta[,type_col] == Disease[1], taxa][1] =1
    }
    if(length(group2) == length(group2[group2<1])) {
      RowFreqMeta[RowFreqMeta[,type_col] == Disease[2], taxa][1] =1
    }

    ####use frequency table to do the feature selection/wilcoxon test
    Genera = RowFreqMeta[,taxa]
    non_zero = sum( Genera!= 0)
    sample_size= length(Genera)
    ratio = non_zero/sample_size  #how many samples have this taxa in the total samples
    #remove low abundance taxa, i.e., if its average relative abundance in either class is <0.5% and its prevelence is at least 20%


    tryCatch({
      fit_glmm1  = wilcox.test(Genera ~ RowFreqMeta[,type_col], data=RowFreqMeta)
      taxalist[[colnames(RowFreqMeta)[taxa]]] <- as.matrix(t(c(fit_glmm1$p.value, MeanType1, SDType1, MeanType2, SDType2)))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
} #end of for loop taxa

  P_value <- data.frame(t(sapply(taxalist,function(x) x)))
  colnames(P_value) = c("P_Wilcoxon", "mean(Type 1)", "std(Type 1)","mean(Type 2)", "std(Type 2)")
  P_value$AdjustP_Wilcoxon = p.adjust(as.numeric(as.character(P_value$P_Wilcoxon)), method = "BH")
  P_value = P_value[as.numeric(as.character(P_value$P_Wilcoxon)) < Cutoff_pvalue,]
##### ratio of mean *(-log10(P))
  meanratio <- vector()
  for (ele in 1:nrow(P_value)) {
    if(P_value$`mean(Type 1)`[ele] > P_value$`mean(Type 2)`[ele]) {
      meanratio[ele] = P_value$`mean(Type 1)`[ele]/P_value$`mean(Type 2)`[ele]
    }else{
      meanratio[ele] = P_value$`mean(Type 2)`[ele]/P_value$`mean(Type 1)`[ele]
    }
    if (meanratio[ele] == Inf) {
      meanratio[ele] <- 200
    }
  }
  NewTable = cbind(P_value, meanratio)
  NewTable2 = cbind(NewTable, NewTable$meanratio*(-log10(as.numeric(as.character(NewTable$P_Wilcoxon)))))
  colnames(NewTable2) = c(colnames(NewTable),"NewScore")


 # return(list(Feature=P_value[order(as.numeric(as.character(P_value$P_Wilcoxon))),], CountData=training))
  return(list(Feature=NewTable2[order(as.numeric(as.character(NewTable2$NewScore)), decreasing = TRUE),], CountData=RowFreqMeta))

}


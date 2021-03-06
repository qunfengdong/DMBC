#' Predict probability from DMBC model
#'
#' @description
#' Based on a set of optimized features from training set, this function predicts the posterior probability for two class labels.
#'
#' @param data Validation dataset with rows are samples, columns are features. The first column should be the sample ID, second column group variable (Disease type, the label you want to classify on).
#' @param testSet Lable unknown testSet without sample IDs.
#' @param auc_out Output object of Cal_AUC() from a validation set.
#' @param col_start An index indicating at which column is the beginning of bacteria (features) data in the validation set. The default is the 3rd column.
#' @param type_col An index indicating at which column is group/type variable in the validation set. The default is the 2nd column.
#' @param Prior1 Prevalence of label1 according to literature or experience. Default is 0.5.
#' @param Prior2 Prevalence of label2 according to literature or experience. Default is 0.5.
#'
#' @return
#'
#' @examples
#' #load the DMBC library
#' library(DMBC)
#'
#' #load training dataset
#' data(training)
#'
#' #load test dataset
#' data(test)
#'
#'## calculate AUC based on training set using 10-fold cv ##
#' auc_out <- Cal_AUC(tfcv(training))
#'
#'## calculate AUC based on training set using leave-one-out cv ##
#' auc_out <- Cal_AUC(loocv(training))
#'
#' #predict unknown test set using training set and auc results.
#' dmbc_predict(data=training,testSet=test,auc_out=auc_out)
#' @export



dmbc_predict <- function(data=data,testSet = testSet,auc_out=auc_out,col_start =3,type_col=2,Prior1=0.5,Prior2=1-Prior1){

  Count= data[,3:ncol(data)]*10000
  Count_Meta = cbind(data[,1:2],Count)
  data = Count_Meta

  testSet = testSet*10000

  namelist <- auc_out$Features[which.max(auc_out$AUC_Area)]
  NameList <- as.vector(unlist(strsplit(as.character(namelist),";")))
  Disease <- levels(data[,type_col])

  NewDF <- as.data.frame(data[,colnames(data) %in% NameList])
  colnames(NewDF) <- colnames(data)[colnames(data) %in% NameList]
  NewDFTotal = data.frame(data[, 1:(col_start-1)],NewDF)

  rep2_Type1 = NewDFTotal[NewDFTotal[,type_col] == Disease[1],]
  rep2_Type2 = NewDFTotal[NewDFTotal[,type_col] == Disease[2],]

  for (taxa in NameList){
    if(nrow(rep2_Type1) == sum(rep2_Type1[,taxa]<1)) {
      rep2_Type1[,taxa][1] =1
    }
    if(nrow(rep2_Type2) == sum(rep2_Type2[,taxa]<1)) {
      rep2_Type2[,taxa][1] =1
    }
  }

  fit3 <- dirmult(round(rep2_Type1[,-(1:(col_start-1))]),epsilon=10^(-4),trace=FALSE)
  fit4 <- dirmult(round(rep2_Type2[,-(1:(col_start-1))]),epsilon=10^(-4),trace=FALSE)

  alpha_Type1 = fit3$gamma
  alpha_Type2 = fit4$gamma

  ###### predict test set #####
  ###### merge in test row

  for (t in NameList){
    if (! (t %in% colnames(testSet)) ){
      testSet <- data.frame(testSet,rep(0,nrow(testSet)))
      colnames(testSet)[ncol(testSet)] <- t
    }
  }

  ###### merge in test row

  NewTestSignature = testSet[, colnames(testSet) %in% NameList]
  NewTestTotal = data.frame(NewTestSignature)

  test_res <- list()
  for (r in 1:nrow(testSet)){
    #print(paste0("testrow: ", r))
    ##### Calculate log of Dirichlet multinomial probability mass function P(x|Type1)
    pdfln_Type1 <- ddirmn(NewTestTotal[r,], t(as.matrix(alpha_Type1)))
    lh_Type1=exp(pdfln_Type1)
    lhP_Type1 = lh_Type1*Prior1

    ##### Calculate log of Dirichlet multinomial probability mass function P(x|Type2)
    pdfln_Type2 <- ddirmn(NewTestTotal[r,], t(as.matrix(alpha_Type2)))
    lh_Type2=exp(pdfln_Type2)
    lhP_Type2 = lh_Type2*Prior2

    Pos_Type1 = lh_Type1/(lh_Type1+lh_Type2)
    Pos_Type2 = lh_Type2/(lh_Type1+lh_Type2)

    PosP_Type1 = lhP_Type1/(lhP_Type1+lhP_Type2)
    PosP_Type2 = lhP_Type2/(lhP_Type1+lhP_Type2)

    test_res[[r]] <- as.matrix(t(c(rownames(testSet)[r],Disease[1], PosP_Type1, Disease[2], PosP_Type2,paste(NameList,collapse=";"))))

  }

  out <- data.frame(t(sapply(test_res,function(x) x)))
  colnames(out) <- c("test_idx","Group1","Group1_prb","Group2","Group2_prb","Features")

  return(out)
}

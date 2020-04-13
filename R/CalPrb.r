#' @title Calculate Likelihood based on Dirichlet-multinomial distribution estimated parameters
#'
#' @description
#' This function estimates parameters from Dirichlet-multinomial distribution.
#' @param FS_out An object from the FeatureSelection()
#' @param testSet A test set in data frame or matrix form. The colnames should have the same bacteria (features) as in the training set.
#' @param col_start An index indicating at which column is the beginning of bacteria (features) data. The default is the 3rd column.
#' @param type_col An index indicating at which column is group/type variable. The default is the 2nd column.
#' @param HighestRank The top number of features inclueded in model. The default is all the features left after filtering.
#' @return A data frame with 17 columns, each row represents a model estimation output.
#' @export
#' @examples
#' data(training)
#'
#' #### Take one row as testSet ####
#' idx <- sample(1:nrow(training),1)
#' test <- training[idx,]
#' train <- training[-idx,]
#'
#' CalPrb(FS(train),test) # This may take up to one minute

CalPrb <- function(FS_out=FS_out,testSet=testSet,col_start=3,type_col=2,HighestRank=nrow(FS_out$Feature)){

  Genus = FS_out$CountData
  SortP = FS_out$Feature

  Disease = levels(Genus[,type_col])

  TotalType1 = nrow(Genus[Genus[,type_col] == Disease[1], ])
  TotalType2 = nrow(Genus[Genus[,type_col] == Disease[2], ])
  Prior1 = TotalType1/(TotalType1+TotalType2)
  Prior2 = TotalType2/(TotalType1+TotalType2)

  lh <- list()
  for(rank in 3:HighestRank) {
   print(rank)

    ######## choose the signature taxa and merge the training data
    ####select the NameList based on the rank
    NameList = as.vector(rownames(SortP)[1:rank])
    NewDF <- Genus[,colnames(Genus) %in% NameList]


    ####In the case of picked 63 bacterium, we don't need to consider others(sum of rest of the columns as others)
    NewDFTotal = data.frame(Genus[, 1:(col_start-1)],NewDF)

    ###### merge in test row

    NewTestSignature = testSet[, colnames(Genus) %in% NameList]
    NewTestTotal = data.frame(testSet[, 1:(col_start-1)],NewTestSignature)

    rep2_Type1 = NewDFTotal[NewDFTotal[,type_col] ==Disease[1],]
    rep2_Type2 = NewDFTotal[NewDFTotal[,type_col] ==Disease[2],]


    ############Estimate Dirchlet-Multinomial parameters
    ###### Estimate the parameters from the control data
    tryCatch({

    fit3 <- dirmult(round(rep2_Type1[,-(1:(col_start-1))]),epsilon=10^(-4),trace=FALSE)
    ###### Estimate the paramenters from the baseline data
    fit4 <- dirmult(round(rep2_Type2[,-(1:(col_start-1))]),epsilon=10^(-4),trace=FALSE)

    #########Calculate the likelihood of the test sample being Type1 and Type2
    alpha_Type1 = fit3$gamma
    alpha_Type2 = fit4$gamma

    testlist <- list()
    for (r in 1:nrow(testSet)){
      ##### Calculate log of Dirichlet multinomial probability mass function P(x|Type1)
      pdfln_Type1 <- ddirmn(round(NewTestTotal[r,-(1:(col_start-1))]), t(as.matrix(alpha_Type1)))
      lh_Type1=exp(pdfln_Type1)
      lhP_Type1 = lh_Type1*Prior1

      ##### Calculate log of Dirichlet multinomial probability mass function P(x|Type2)
      pdfln_Type2 <- ddirmn(round(NewTestTotal[r,-(1:(col_start-1))]), t(as.matrix(alpha_Type2)))
      lh_Type2=exp(pdfln_Type2)
      lhP_Type2 = lh_Type2*Prior2

      Pos_Type1 = lh_Type1/(lh_Type1+lh_Type2)
      Pos_Type2 = lh_Type2/(lh_Type1+lh_Type2)

      PosP_Type1 = lhP_Type1/(lhP_Type1+lhP_Type2)
      PosP_Type2 = lhP_Type2/(lhP_Type1+lhP_Type2)

      #### Create truth labels ####

      if(testSet[r,type_col] == Disease[1]) {
        Type1Label = 1
        Type2Label = 0
      } else if (testSet[r,type_col] == Disease[2]) {
        Type1Label = 0
        Type2Label = 1
      }

      testlist[[r]] <- as.matrix(t(c(Disease[1], Disease[2], rownames(testSet)[r], rank,  PosP_Type1, PosP_Type2, Type1Label, Type2Label,paste(NameList,collapse=";"))))
    }
    lh[[rank]] <- data.frame(t(sapply(testlist,'[')))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  } #end of rank

  lh_table <- data.frame(do.call(rbind,lh))
  colnames(lh_table) = c("Type1", "Type2", "row", "feature_rank", "Type1_Posterior_Prb","Type2_Postereior_Prb", "Type1Label", "Type2Label","SelectedFeatures" )

  return(lh_table)
}

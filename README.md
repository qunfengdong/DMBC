Dirichlet-Multinomial Bayes Classifier (DMBC) for Microbiome Classification
---------------------------------------------------------------------------
This package implements the machine-learning method described by Gao et al (2016) for microbiome classification using a Bayes classifier based on the Dirichlet-Multinomial distribution.  In addition to classification, the package also identifies a subset of microbial taxa that can achieve the maximum classification accuracy.

## Citation
Xiang Gao, Huaiying Lin, Qunfeng Dong (2017); [A Dirichlet-Multinomial Bayes Classifier for Disease Diagnosis with Microbial Compositions](http://msphere.asm.org/content/msph/2/6/e00536-17.full.pdf), mSphere, Volume: 2, Issue: 6.

## Update
* 12/11/2018	v1.1.0	changes include 1) Relative abaundance table is used as input. 2) Variable "Others", which is the sum of all the variables except for the selected features, is not considered as a feature in this version. 3) Ten fold validation function is temporarily disabled. 4) Pi score instead of the wilcoxon p value is used to rank the importantce of features (Xiao Y, Hsiao TH, Suresh U, et al. A novel significance score for gene selection and ranking. Bioinformatics. 2012;30(6):801-7. )
* 3/12/2018 v0.2.1 Increase consistency of feature selection.
* 3/7/2018 v0.2.0 Add ten-fold cross-validation function in addition to the existing leave-one-out cross-validation function

## Installation

Please make sure you have devtools installed in your R and do the following:

```R
library(devtools)
install_github("qunfengdong/DMBC")
```

## Usage

There are totally 6 functions included in the package. The most important one is dmbc_predict, which will predict probability for a test set from DMBC model given a training set.

* best_cm	: Confusion Matrix for the Best Model
* CalPrb	: Calculate Likelihood based on Dirichlet-multinomial distribution estimated parameters
* Cal_AUC	: Calculate Area Under the ROC Curve from Leave-One-Out Cross-Validation
* dmbc_predict	: Predict probability from DMBC model 
* FS	: Filter and Feature Selection based on Wilcoxon Rank-Sum test
* loocv	: Leave-one-out cross-validation
* test	: Example Testing Data
* training	: Example Training Data

## Contributors

* Xiang Gao (Author)
* Qunfeng Dong (Author)
* Huaiying Lin (Maintainer)

## License

MIT

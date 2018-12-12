Dirichlet-Multinomial Bayes Classifier (DMBC) for Microbiome Classification
---------------------------------------------------------------------------
This package implements the machine-learning method described by Gao et al (2016) for microbiome classification using a Bayes classifier based on the Dirichlet-Multinomial distribution.  In addition to classification, the package also identifies a subset of microbial taxa that can achieve the maximum classification accuracy.

## Citation
Xiang Gao, Huaiying Lin, Qunfeng Dong (2017); [A Dirichlet-Multinomial Bayes Classifier for Disease Diagnosis with Microbial Compositions](http://msphere.asm.org/content/msph/2/6/e00536-17.full.pdf), mSphere, Volume: 2, Issue: 6.

## Update
* 3/7/2018 v0.2.0 Add ten-fold cross-validation function in addition to the existing leave-one-out cross-validation function
* 3/12/2018 v0.2.1 Increase consistency of feature selection.
* 12/11/2018	v1.1.0	changes include 1) use relative abaundance as input, and in the function, it will be multiple by 10000 and rounding to count data. 2) Will not sum the rest of the variable to create a new vailable called "rest". 3) use Pi score instead of the wilcoxon p value to rank the importantce of features (Xiao Y, Hsiao TH, Suresh U, et al. A novel significance score for gene selection and ranking. Bioinformatics. 2012;30(6):801-7. )

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
* tfcv	: Ten-fold cross-validation
* test	: Example Testing Data
* training	: Example Training Data

## Contributors

* Xiang Gao (Author)
* Qunfeng Dong (Author)
* Huaiying Lin (Maintainer)

## License

MIT

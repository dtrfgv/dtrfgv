
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tplda

Attempt to create an R package from prototype available at
[“Classification tree algorithm for grouped
variables”](https://github.com/apoterie/TPLDA).

## Installation

Don’t try to install this package :sweat\_smile:

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(tplda)

data(data_pour_tester_fonctions)

train<-data[which(data[,1]=="train"),-1]           # negative index into the `data` 
test<-data[which(data[,1]=="test"),-1]             # object specifying all rows and all columns 
validation<-data[which(data[,1]=="validation"),-1] # except the first column.

forest<-rfgv(train,
             group=group,
             groupImp=group,
             ntree=1,
             mtry_group=3,
             sampvar=TRUE,
             sampvar_type=2,
             maxdepth=2,
             kfold=3,
             replace=T,
             case_min=1,
             sampsize=nrow(train),
             mtry_var=rep(2,5),
             grp.importance=TRUE,
             test=test,
             keep_forest=F,
             crit=1,
             penalty="No")

print(forest$importance)
#>    MeanDecrAcc MeanDecrAccNor
#> 1 -0.008474576   -0.001694915
#> 2  0.084745763    0.016949153
#> 3 -0.042372881   -0.008474576
#> 4 -0.016949153   -0.003389831
#> 5  0.042372881    0.008474576
```

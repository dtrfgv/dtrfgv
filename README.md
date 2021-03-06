
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Decision Trees and Random Forests for Grouped Variables

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/dtrfgv)](https://cran.r-project.org/package=dtrfgv)
[![Build
Status](https://travis-ci.org/dtrfgv/dtrfgv.svg?branch=master)](https://travis-ci.org/dtrfgv/dtrfgv)
[![Coverage
status](https://codecov.io/gh/dtrfgv/dtrfgv/branch/master/graph/badge.svg)](https://codecov.io/github/dtrfgv/dtrfgv?branch=master)
<!-- badges: end -->

Attempt to create an R package from prototype available at
[“Classification tree algorithm for grouped
variables”](https://github.com/apoterie/TPLDA).

## Installation

Don’t try to install this package :sweat\_smile:

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(dtrfgv)

data(rfgv_dataset)
data(group)

data  <- rfgv_dataset 
train <- data[which(data[,1]=="train"),-1]           # negative index into the `data` 
test  <- data[which(data[,1]=="test"),-1]             # object specifying all rows and all columns 
validation<-data[which(data[,1]=="validation"),-1] # except the first column.

forest<-rfgv(train,
             group=group,
             groupImp=group,
             ntree=4,
             mtry_group=3,
             sampvar=TRUE,
             maxdepth=2,
             replace=TRUE,
             case_min=1,
             sampsize=nrow(train),
             mtry_var=rep(2,5),
             grp.importance=TRUE,
             test=test,
             keep_forest=FALSE,
             crit=1,
             penalty="No")
  
print(forest$importance)
#>   MeanDecrAcc MeanDecrAccNor
#> 1  0.01260190    0.002520380
#> 2  0.19578711    0.039157422
#> 3 -0.01797763   -0.003595526
#> 4 -0.01288626   -0.002577253
#> 5 -0.01007369   -0.002014739
```

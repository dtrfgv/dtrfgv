
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Tree Penalized Linear Discriminant Analysis

[![Build
Status](https://travis-ci.org/pnavaro/tplda.svg?branch=master)](https://travis-ci.org/pnavaro/tplda)
[![Coverage
status](https://codecov.io/gh/pnavaro/tplda/branch/master/graph/badge.svg)](https://codecov.io/github/pnavaro/tplda?branch=master)

Attempt to create an R package from prototype available at
[“Classification tree algorithm for grouped
variables”](https://github.com/apoterie/TPLDA).

## Installation

Don’t try to install this package :sweat\_smile:

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(tplda)

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
#> 1  0.04328686    0.008657372
#> 2  0.19115480    0.038230960
#> 3  0.01059353    0.002118706
#> 4  0.02221928    0.004443857
#> 5  0.02163947    0.004327894
```

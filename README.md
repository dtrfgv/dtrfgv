
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
#> 1 -0.007518797   -0.001503759
#> 2  0.105263158    0.021052632
#> 3  0.000000000    0.000000000
#> 4 -0.075187970   -0.015037594
#> 5  0.000000000    0.000000000
```

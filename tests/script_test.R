rm(list=ls())

library(tplda)
library(gtools)
library(rpart)
library(e1071)
library(fields)
library(MASS)
library(UBL)
library(ElemStatLearn)
library(rpart.plot)
library(randomForest)
library(pryr)
require(compiler)

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

forest$importance
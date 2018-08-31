rm(list=ls())

chemin<-'/Users/navaro/R-projects/phoneme/tests/Test_Pierre/'
## CHARGEMENT DES FONCTIONS
load(paste(chemin,'data_pour_tester_fonctions.RData',sep=""))

source(paste(chemin, 'Fonctions_RF_evo_Modif_pierre.R', sep = ""))
source(paste(chemin, 'Fonctions_performance.R', sep = ""))
source(paste(chemin, 'Fonctions_CARTGVRF_sampvar1.R', sep = ""))# CARTGV original et variante avec sampvar_type=1
source(paste(chemin, 'Fonctions_CARTGVRF_sampvar2.R', sep = ""))# variante de CARTGV avec sampvar_type=2
source(paste(chemin, 'Fonctions_RPART_sup.R', sep = ""))

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

## EXECUTABLE


train<-data[which(data[,1]=="train"),-1]
test<-data[which(data[,1]=="test"),-1]
validation<-data[which(data[,1]=="validation"),-1]


forest<-rfgv(train,
             group=group,
             groupImp=group,
             ntree=10,
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



# rm(list=ls())
# 
# library(here)
# 
# load(here('data','data_pour_tester_fonctions.RData'))
# 
# source(here('R', 'Fonctions_RF_evo_Modif_pierre.R'))
# source(here('R', 'Fonctions_performance.R'))
# source(here('R', 'Fonctions_CARTGVRF_sampvar1.R'))# CARTGV original et variante avec sampvar_type=1
# source(here('R', 'Fonctions_CARTGVRF_sampvar2.R'))# variante de CARTGV avec sampvar_type=2
# source(here('R', 'Fonctions_RPART_sup.R'))
# 
# library(gtools)
# library(rpart)
# library(e1071)
# library(fields)
# library(MASS)
# library(UBL)
# library(ElemStatLearn)
# library(rpart.plot)
# library(randomForest)
# library(pryr)
# require(compiler)
# 
# train<-data[which(data[,1]=="train"),-1]
# test<-data[which(data[,1]=="test"),-1]
# validation<-data[which(data[,1]=="validation"),-1]
# 
# 
# forest<-rfgv(train,
#              group=group,
#              groupImp=group,
#              ntree=10,
#              mtry_group=3,
#              sampvar=TRUE,
#              sampvar_type=2,
#              maxdepth=2,
#              kfold=3,
#              replace=T,
#              case_min=1,
#              sampsize=nrow(train),
#              mtry_var=rep(2,5),
#              grp.importance=TRUE,
#              test=test,
#              keep_forest=F,
#              crit=1,
#              penalty="No")
# 
# 
# forest$importance
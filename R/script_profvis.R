rm(list=ls())

library(tplda)
library(profvis)

profvis({
data(data_pour_tester_fonctions)

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
             replace=TRUE,
             case_min=1,
             sampsize=nrow(train),
             mtry_var=rep(2,5),
             grp.importance=TRUE,
             test=test,
             keep_forest=FALSE,
             crit=1,
             penalty="No")
})

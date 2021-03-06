context("rfgv")

test_that("rfgv with a random forest of 4 trees ", {
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
  
  
  expect_equal(which.max(forest$importance[[1]]), 2)
})

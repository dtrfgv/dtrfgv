############################################################################################
#################    Random Forests algorithm for grouped variables     ####################
############################################################################################


#' entropy
#' 
#' Calculate the value of the entropy for binary classification (i.e. Y= 0 or 1)
#' 
#' @param p a vector containit the probabilities '\code{P[Y=0]}' and '\code{P[Y=1]}'
#'
#' @return the value of the entropy
#'
entropy <- function(p) {
  if (any(p == 1)) return(0)# -sum(p *
  -sum(p*log(p,2)) 
}

#' gini
#'
#' Calculate the value of the Gini index for binary classification (i.e. Y= 0 or 1)
#'
#' @param p p a vector containit the probabilities '\code{P[Y=0]}' and '\code{P[Y=1]}'
#'
#' @return the value of the Gini index
#'
gini <- function(p) {
  sum(p * (1 - p))
}

#' bsamples
#' 
#' Generates the bootstrap samples used to built a RFGV forest
#' 
#'
#' @param ntree number of boostrap samples (i.e. the number of trees in the RFGV forest)
#' @param data  the learning data set. it must be a data frame
#' @param sampsize an interger indicating the size of the boostrap samples
#' @param replace a boolean indicating if sampling of cases is done with or without replacement?
#' 
#' @return a list with elements
#'          - bsamples:  a matrix which with '\code{sampsize}' lines and '\code{B}' columns which contains the indices of the observations
#'                       belonging to each boostrap sample
#'          - oobsamples: a list that contains the indices of the observations belonging to each out-of-bag sample
#'
bsamples<-function(ntree,data,sampsize,replace){
  N<-nrow(data)
  bsamples<-apply(matrix(rep(1:N,ntree),ncol=ntree,nrow=N),
                  2,sample,size=sampsize,replace=replace)
  oobsamples<-lapply(1:ntree,function(x)setdiff(1:N,bsamples[,x]))
  return(list(bsamples=bsamples,
              oobsamples=oobsamples))
}


#' group.selection
#' 
#' Selects randomly '\code{mtry}' groups
#' 
#'
#' @param group a vector indicating the group label of each variable
#' @param mtry  a number indicating the number of selected group
#'
#' @return groups a vector indicating the label of the selected groups
#'
group.selection<-function(group,
                          mtry=sqrt(unique(group[!is.na(group)])))
  {
  
  if(length(which(is.na(group)))>0){
    print("Warning: there are NA in group; The fucntion will ignore these NA")
    label<-unique(group[!is.na(group)])
  }else{
    label<-unique(group)
  }
  
  p<-length(label)
  groups<-NULL
  
  if(mtry<=p){
    groups<-sample(label,mtry,replace=FALSE)
  }else{
    print("WARNING: mtry must be lower or equal to p")
    groups<-label
  }
  return(groups)
}


#' rfgv
#' 
#' Random Forest for Grouped Variables. Only implement for binary classification. 
#' The function builts a large number of random decision trees based on a variant of the CARTGV method.
#'
#' @param data a data frame containing the response value (for the first variable)  and the predictors 
#' and used to grow the tree. The name of the response value must be "Y". 
#' The response variable must be the first variable of the data frame and the variable meust be coded as the two levels "0" and "1".
#' @param group a vector with the group number of each variable. 
#  (WARNING : if there are "\code{p}" goups, the groups must be numbers from "\code{1}" to "\code{p}" in increasing order. The group label of the response variable is missing (i.e. NA))
#' @param groupImp a vector which indicates the group number of each variable (for the groups used to compute the group importance).
#' @param ntree an integer indicating the number of trees to grow 
#' @param mtry_group an integer the number of variables randomly samples as candidates at each split.
#' @param maxdepth an integer indicating the maximal depth for a split-tree. The default value is 2.
#' @param replace a boolean indicating if sampling of cases is done with or without replacement?
#' @param sampsize an interger indicating the size of the boostrap samples.
#' @param case_min an integer indicating the minimun number of cases/non cases in a terminal nodes. The default is 1.
#' @param grp.importance a boolean indicating if the importance of each group need to be computed
#' @param test an independent data frame containing the same variables that "\code{data}".
#' @param keep_forest a boolean indicating if the forest will be retained in the output object
#' @param crit an integer indicating the impurity function used (1=Gini index / 2=Entropie/ 3=Misclassification rate)
#' @param penalty a boolean indicating if the decrease in node impurity must take account of the group size. Four penalty are available: "No","Size","Root.size" or "Log".
#' @param sampvar a boolean indicating if within each splitting tree, a subset of variables is drawn for each group
#' @param mtry_var a vector of length the number of groups. It indicates the number of drawn variables for each group. Usefull only if sampvar=TRUE 
#'
#' @return a list with elements:
#'         - predicted: the predicted values of the observations in the training set named "\code{data}". The i-th element being the prediction from the ith tree and based
#'                      on the i-th out-of-bag sample. The i-th element is missing if the i-th observation is not part of the the i-th out-of-bag sample.
#'         - importance: a data frame with two coloums. The first column provides the value of the permutation importance of each group 
#'                       and the second one gives the value of the permutation importance of each group normalized by the size of the group
#'         - err.rate:   a vector error rates of the prediction on the training set named "\code{data}", the i-th element being the (OOB) error rate 
#'                       for all trees up to the i-th.            
#'         - vote: a data frame with one row for each input data point and one column for each class ("0" and "1", in this order), giving the fraction 
#'                 number of (OOB) ‘votes’ from the random forest.
#'         - pred: the predicted values of the observations in the training set named "\code{data}". It correspond to the majority vote computed by using the 
#'                 matrix of predictions "\code{predicted}".
#'         - confusion:  the object returned by the function "\code{xtab_function}". There are the confusion matrix of the prediction (based on OOB data) and
#'                       the associated statistics. For more details, see the function "\code{xtab_function}".
#'         - err.rate.test:  (Only if test!=NULL) a vector error rates of the prediction on the test set named "\code{test}", the i-th element being the error rate 
#'                           for all trees up to the i-th. 
#'         - vote.test:  (Only if test!=NULL) a data frame with one row for each observtion in "\code{test}" and one column for each class ("0" and "1", in this order), 
#'                       giving the number of ‘votes’ from the random forest.
#'         - pred.test:  (Only if test!=NULL) the predicted values of the observations in "\code{test}".
#'         - confusion.test:  (Only if test!=NULL)  the object returned by the function "\code{xtab_function}". There are the confusion matrix of the prediction (based on 
#'                            "\code{test}") and the associated statistics. For more details, see the function "\code{xtab_function}".
#'         - oob.times:  number of times that an observation in the training set named "\code{data}" is ‘out-of-bag’ (and thus used in computing OOB error estimate)
#'         - keep_forest:  a boolean indicating if the forest will be retained in the output object
#'         - sampvar:  a boolean indicating if within each splitting tree, a subset of variables is drawn for each group
#'         - maxdepth:  an integer indicating the maximal depth for a split-tree. The default value is 2
#'         - mtry_group:  an integer the number of variables randomly samples as candidates at each split
#'         - mtry_var: a vector of length the number of groups. It indicates the number of drawn variables for each group. Usefull only if sampvar=TRUE
#'         - ntree: an integer indicating the number of trees to grow.
#'                       
#' @export
#'
#' @examples
#' data(rfgv_dataset)
#' data(group)
#' data <- rfgv_dataset
#' train<-data[which(data[,1]=="train"),-1]           # negative index into the `data` 
#' test<-data[which(data[,1]=="test"),-1]             # object specifying all rows and all columns 
#' validation<-data[which(data[,1]=="validation"),-1] # except the first column.
#' 
#' forest<-rfgv(train,
#'              group=group,
#'              groupImp=group,
#'              ntree=1,
#'              mtry_group=3,
#'              sampvar=TRUE,
#'              maxdepth=2,
#'              replace=TRUE,
#'              case_min=1,
#'              sampsize=nrow(train),
#'              mtry_var=rep(2,5),
#'              grp.importance=TRUE,
#'              test=test,
#'              keep_forest=FALSE,
#'              crit=1,
#'              penalty="No")
rfgv<-function(data,
               group,
               groupImp,
               ntree=200,
               mtry_group=floor(sqrt(length(unique(group[!is.na(group)])))),
               maxdepth=1,
               replace=T,
               sampsize=ifelse(replace==T,nrow(data),floor(0.632*nrow(data))),
               case_min=1,
               grp.importance=TRUE,
               test=NULL,
               keep_forest=F,
               crit=1,
               penalty="No",
               sampvar=FALSE,
               mtry_var)
{
  nbgrp<-length(unique(group[!is.na(group)]))
  predicted<-matrix(rep(NA,nrow(data)*ntree),ncol=ntree)
  vote<-matrix(rep(NA,nrow(data)*2),ncol=2)
  err<-rep(0,ntree)
  if(!is.null(test)){
    vote.test<-matrix(rep(NA,nrow(test)*2),ncol=2)
    test_predicted<-matrix(rep(NA,nrow(test)*ntree),ncol=ntree)
    err.test<-rep(0,ntree)
  }
  if(keep_forest==TRUE){
    forest<-list()
  }else{
    forest<-NULL
  }
  if(grp.importance==TRUE){
    nbgrpImp<-length(unique(groupImp[!is.na(groupImp)]))
    acc<-rep(NA,1)
    decreaacc<-matrix(rep(0,ntree*length(unique(groupImp[!is.na(groupImp)]))),nrow=ntree)
    MeanDecrAcc<-rep(0,nbgrpImp)
  }
  samples<-bsamples(ntree=ntree,data,sampsize,replace)
  for(b in 1:ntree){
      tree<-cartgv.rf(data[unlist(samples$bsamples[,b]),],
                      group,
                      crit=crit,
                      case_min=case_min,
                      maxdepth=maxdepth,
                      p=mtry_group,
                      penalty=penalty,
                      mtry_var=mtry_var)
      
      predicted[unlist(samples$oobsamples[[b]]),b]<-as.numeric.factor(predict_cartgv.rf(data[unlist(samples$oobsamples[[b]]),], 
                                                                                        tree$tree,
                                                                                        tree$tree_split)$hat.Y)
      
    err[b]<-length(which(predicted[unlist(samples$oobsamples[[b]]),b]!=data$Y[unlist(samples$oobsamples[[b]])]))/length(data$Y[unlist(samples$oobsamples[[b]])])
    if(grp.importance==TRUE){
      acc<-1-as.numeric(impurity.cartgv.rf(data[unlist(samples$oobsamples[[b]]),],
                                             list(tree$tree),tree)$impurete$Misclass[1])
      groupselec<-unique(as.numeric.factor(tree$tree$var[which(tree$tree$action == "1")]))
      groupselecImp<-unique(groupImp[which(group%in%groupselec)])
      decreaacc[b,groupselecImp]<-as.numeric(unlist(sapply(groupselecImp,
                                                          function(x)grpimpperm.rf(num_group=x,
                                                                                   data=data,
                                                                                   oobsamples=samples$oobsamples[[b]],
                                                                                   group=groupImp,
                                                                                   tree=tree,
                                                                                   impurityacc=acc))))
      acc<-NULL
    }
    if(!is.null(test)){
        test_predicted[,b]<-as.numeric.factor(predict_cartgv.rf(test,
                                                                tree$tree,
                                                                tree$tree_split)$hat.Y)
      
         err.test[b]<-length(which(test_predicted[,b]!=test$Y))/length(test$Y)
      
    }
    if(keep_forest==TRUE){
      forest[[b]]<-tree
    }
    rm(tree)
    gc()
  }
  err.rate<-apply(cbind(cumsum(err),1:ntree),1,function(x)x[1]/x[2])
  vote[,2]<-rowMeans(predicted,na.rm=T)
  vote[,2]<-ifelse(!is.nan(vote[,2]),vote[,2],NA)
  vote[,1]<-ifelse(!is.na(vote[,2]),1-vote[,2],NA)
  vote<-as.data.frame(vote)
  names(vote)<-c("hat.Y=0","hat.Y=1")
  pred<-ifelse(vote[,1]>vote[,2],"0",ifelse(vote[,1]<=vote[,2],"1",NA))
  oob.times<-rowSums(ifelse(!is.na(predicted),1,0),na.rm=T)
  confusion<-xtab_function(pred[!is.na(pred)],data$Y[!is.na(pred)])
  
  if(!is.null(test)){
    err.rate.test<-apply(cbind(cumsum(err.test),1:ntree),1,function(x)x[1]/x[2])
    vote.test[,2]<-rowMeans(test_predicted,na.rm=T)
    vote.test[,1]<-ifelse(!is.na(vote.test[,2]),1-vote.test[,2],NA)
    vote.test<-as.data.frame(vote.test)
    names(vote.test)<-c("hat.Y=0","hat.Y=1")
    pred.test<-ifelse(vote.test[,1]>vote.test[,2],"0",ifelse(vote.test[,1]<=vote.test[,2],"1",NA))
    confusion.test<-xtab_function(pred.test,test$Y)
  }else{
    err.rate.test<-NULL
    pred.test<-NULL
    confusion.test<-NULL
    vote.test<-NULL
  }
  if(grp.importance==TRUE){
    MeanDecrAcc<-colSums(decreaacc)/ntree
    MeanDecrAccNor<-apply(cbind(MeanDecrAcc,
                                as.numeric(table(groupImp[-1]))),
                          1,
                          function(x)x[1]/x[2])
    importance<-as.data.frame(cbind(MeanDecrAcc,MeanDecrAccNor))
  }else{
    importance<-NULL
  }
  
  return(list(forest=forest,
              predicted=predicted,
              importance=importance,
              err.rate=err.rate,
              vote=vote,
              pred=pred,
              confusion=confusion,
              err.rate.test=err.rate.test,
              vote.test=vote.test,
              pred.test=pred.test,
              confusion.test=confusion.test,
              oob.times=oob.times,
              keep_forest=keep_forest,
              sampvar=sampvar,
              maxdepth=maxdepth,
              mtry_group=mtry_group,
              mtry_var=mtry_var,
              ntree=ntree))
}



#' surrogate_split
#'
#' @param Ystar
#' @param node
#' @param group
#' @param igroup
#' @param penalty
#' @param label
#' @param maxdepth
#'
#' @return
#' @export
#'
surrogate_split<-function(Ystar,node,group,igroup,penalty='No',label,maxdepth){
  ivar<-which(group==igroup)
  group.size<-length(ivar)
  data<-as.data.frame(cbind(Ystar,node[,ivar]))
  names(data)<-c("Ystar",names(node)[ivar])
  Ypred<-rep(NA,length(Ystar))
  Gain_Gini<-0
  Gain_Ent<-0
  Gain_Clas<-0
  formula<-as.formula(paste("Ystar ~ ",paste(names(data)[-1],collapse="+")))
  cart <- rpart(formula,data,control=rpart.control(minsplit=1,cp=-0.999,maxdepth = maxdepth,maxsurrogate =0,maxcompete = 0 ),
                method = "class",parms = list(split = "gini"))
  Ystarpred<-as.numeric(as.character(predict(cart, type = "class", newdata=data)))
  modalities<-unique(as.numeric(as.character(cart$where)))
  for(k in 1:length(modalities)){
    Ypred[which(unlist(cart$where)==modalities[k])]<-ifelse(sum(as.numeric(as.character(node$Y[which(unlist(cart$where)==modalities[k])])))/(length(node$Y[which(unlist(cart$where)==modalities[k])]))<0.5,"0","1")
  }
  #if(sum(Ypred==2,na.rm=T)>0){Ypred<-ifelse(Ypred==2,"1","0")}
  tab<-table(Ypred,node$Y)
  Gain_Gini<-length(node$Y)*(gini(prop.table(table(node$Y))))
  Gain_Ent<-length(node$Y)*(entropy(prop.table(table(node$Y))))
  Gain_Clas<-length(node$Y)*((length(which(label!=node$Y)))- length(which(Ypred!=node$Y)))
  for(k in 1:length(modalities)){
    node1<-node[which(unlist(cart$where)==modalities[k]),]
    tab1<-table(Ypred[which(unlist(cart$where)==modalities[k])],node$Y[which(unlist(cart$where)==modalities[k])])
    prop_n1<-dim(node1)[1]/length(node$Y)
    Gain_Gini<-Gain_Gini- length(node$Y)*(prop_n1*gini(prop.table(tab1[1,])))
    Gain_Ent<-Gain_Ent - length(node$Y)*(prop_n1*entropy(prop.table(tab1[1,])))
  }
  if(penalty=="Size"){
    Gain_Gini<-Gain_Gini/group.size
    Gain_Ent<-Gain_Ent/group.size
    Gain_Clas<-Gain_Clas/group.size
  }
  if(penalty=="Root.size"){
    Gain_Gini<-Gain_Gini/sqrt(group.size)
    Gain_Ent<-Gain_Ent/sqrt(group.size)
    Gain_Clas<-Gain_Clas/sqrt(group.size)
  }
  if(penalty=="Log"){
    if(group.size>1){
      Gain_Gini<-Gain_Gini/(log(group.size))
      Gain_Ent<-Gain_Ent/(log(group.size))
      Gain_Clas<-Gain_Clas/log(group.size)
    }
  }
   return(list(Gain_Gini=Gain_Gini,Gain_Ent=Gain_Ent,Gain_Clas=Gain_Clas,Ypred=Ypred))
  }


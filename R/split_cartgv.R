# =========================================================================================
# split.cartgv()
# =========================================================================================


#' Title
#'
#' @param node 
#' @param group 
#' @param igroups 
#' @param label 
#' @param maxdepth 
#' @param penalty 
#' @param sampvar 
#' @param mtry_var 
#'
#' @return
#' @export
#' @importFrom stats predict as.formula
#' @examples
split.cartgv<-function(node,group,igroups,label,maxdepth=2,penalty="No",sampvar="FALSE",
                       mtry_var=sapply(as.numeric(table(group[!is.na(group)])),function(x)floor(sqrt(x)))){
  nb_group<-length(unique(group[which(group%in%igroups)]))
  Gain_Gini<-rep(0,nb_group)
  Gain_Ent<-rep(0,nb_group)
  Gain_Clas<-rep(0,nb_group)
  carts <- list()
  nb_nodes<-rep(0,nb_group)
  pred<-matrix(rep(NA,dim(node)[1]*nb_group),ncol=nb_group)
  var_selec<-rep(0,length(group[!is.na(group)]))
  for(j in 1:nb_group){
    if(sampvar=="TRUE"){
      ivar<-sample(which(group==igroups[j]),size=min(mtry_var[igroups[j]], length(which(group==igroups[j]))),replace=FALSE)
      var_selec[ivar]<-1
    }else{
      ivar<-which(group==igroups[j])
    }
    formula<-as.formula(paste("Y ~ ",paste(names(node)[ivar],collapse="+")))
    cart <- rpart(formula,node,control=rpart.control(minsplit=1,cp=-0.999,maxdepth = maxdepth,maxsurrogate =0,maxcompete = 0 ),method = "class",parms = list(split = "gini"))
    pred[,j]<-as.numeric(as.character(predict(cart, type = "class", newdata=node)))
    carts[[j]] <- cart
    cart<-NULL
  }
  for(j in 1:nb_group){
    if(sum(pred[,j]==2,na.rm=T)>0){pred[,j]<-ifelse(pred[,j]==2,"1","0")}
    tab<-table(pred[,j],node$Y)
    modalities<-unique(carts[[j]]$where)
    Gain_Gini[j]<-length(node$Y)*(gini(prop.table(table(node$Y))))
    Gain_Ent[j]<-length(node$Y)*(entropy(prop.table(table(node$Y))))
    Gain_Clas[j]<-length(node$Y)*((length(which(label!=node$Y)))- length(which(pred[,j]!=node$Y)))
    for(k in 1:length(modalities)){
      node1<-node[which(unlist(carts[[j]]$where)==modalities[k]),]
      tab1<-table(pred[which(unlist(carts[[j]]$where)==modalities[k]),j],node$Y[which(unlist(carts[[j]]$where)==modalities[k])])
      prop_n1<-dim(node1)[1]/length(node$Y)
      Gain_Gini[j]<-Gain_Gini[j]- length(node$Y)*(prop_n1*gini(prop.table(tab1[1,]))) 
      Gain_Ent[j]<-Gain_Ent[j] - length(node$Y)*(prop_n1*entropy(prop.table(tab1[1,])))
    }
    if(penalty=="Size"){
      group.size<-length(which(group==igroups[j]))
      Gain_Gini[j]<-Gain_Gini[j]/group.size
      Gain_Ent[j]<-Gain_Ent[j]/group.size
      Gain_Clas[j]<-Gain_Clas[j]/group.size
    }
    if(penalty=="Root.size"){
      group.size<-length(which(group==igroups[j]))
      Gain_Gini[j]<-Gain_Gini[j]/sqrt(group.size)
      Gain_Ent[j]<-Gain_Ent[j]/sqrt(group.size)
      Gain_Clas[j]<-Gain_Clas[j]/sqrt(group.size)
    }
    if(penalty=="Log"){
      group.size<-length(which(group==igroups[j]))
      if(group.size>1){
        Gain_Gini[j]<-Gain_Gini[j]/(log(group.size))
        Gain_Ent[j]<-Gain_Ent[j]/(log(group.size))
        Gain_Clas[j]<-Gain_Clas[j]/log(group.size)
      }
      
    }
  }
  return(list(Gain_Gini=Gain_Gini, Gain_Ent=Gain_Ent, Gain_Clas=Gain_Clas, carts=carts, pred=pred,igroups=igroups,var_selec=var_selec))
}

############################################################################################
########################             CARTGV original                #######################
########################                   et                       #######################
########################      Variante avec sampvar_type=1          #######################
############################################################################################

library(e1071)### l'indice de rand utilise la fonction classagreement() du package



# =========================================================================================
# perm()
# =========================================================================================

#' Title
#'
#' @param oobsamples 
#' @param data 
#' @param num_group 
#' @param group 
#'
#' @return
#' @export
#' @importFrom gtools permute
#' @examples
perm<-function(oobsamples,data,num_group,group){
  
  data_perm<-data[oobsamples,]
  iperm<-permute(1:nrow(data_perm))
  data_perm[,which(group==num_group)]<-data_perm[iperm,which(group==num_group)]
  names(data_perm)<-names(data)
  
  return(data_perm)
}

# =========================================================================================
# grpimpperm()
# =========================================================================================



#' Title
#'
#' @param num_group 
#' @param data 
#' @param oobsamples 
#' @param group 
#' @param tree 
#' @param impurityacc 
#'
#' @return
#' @export
#'
#' @examples
grpimpperm<-function(num_group,data,oobsamples,group,tree,impurityacc){
  data_perm<-perm(oobsamples,data,num_group,group)
  accperm<-1-impurity.cartgv(data_perm,list(tree$tree),tree)$impurete$Misclass# accurancy=1-misclass
  DecreaseAcc<-impurityacc-accperm
  return(DecreaseAcc)
}


#' Title grpimpgini
#'
#' @param num_group 
#' @param groupselec 
#' @param tree 
#'
#' @return
#' @export
#'
#' @examples
grpimpgini<-function(num_group,groupselec,tree){
  DecreaseImpurity<-0
  times<-which(tree$tree$var==num_group)
  if(num_group%in%groupselec){
    for(i in times ){
      p<-(unlist(as.numeric(as.character(tree$tree$n[i])))
          /unlist(as.numeric(as.character(tree$tree$n[1]))))
      DecreaseImpurity<-DecreaseImpurity
      +(p*unlist(tree$split[[i]]$Gain_Gini[which(tree$split[[i]]$igroups==num_group)]))
    }
  }
  return(DecreaseImpurity)
}



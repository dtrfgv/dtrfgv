############################################################################################
########################           Variante de CARTGV              ########################
########################           avec sampvar_type=2             ########################
############################################################################################


### 
#' Classification Random Tree for Grouped Variables using Random Forest.
#'
#' A utiliser quand sampvar=TRUE and sampvar_type="2 ==> permet de construire un arbre 
#' cartgv tel que dans chaque arbre de 
#' coupure un sous-ensemble de variables est tiré au hasard avant chaque coupure
#' 
#' @param data 
#' @param group 
#' @param crit 
#' @param case_min 
#' @param maxdepth 
#' @param p 
#' @param penalty 
#' @param mtry_var 
#'
#' @return
#' @export
#'
#' @examples
cartgv.rf<-function(data,
                    group,
                    crit=1,
                    case_min=1,
                    maxdepth=2,
                    p=floor(sqrt(length(unique(group[!is.na(group)])))),
                    penalty="No",
                    mtry_var=sapply(as.numeric(table(group[!is.na(group)])),function(x)floor(sqrt(x))))
  {
  ##Initialisation
  tree<-NULL
  parent<-c()#Noeud ancetre du noeud
  depth<-c()#Profondeur du noeuds
  var<-c()#Groupe selectionné pour couper
  carts<-list()#CART définissant la coupure
  n<-c()# Effectif dans le noeud considéré
  n_case<-c()# Effectif de Y=1 dans le noeud considéré
  n_noncase<-c()# Effectif de Y=0 dans le noeud considéré
  yval<-c()#Prédiction de CART = règle de prediction basée sur la classe majoritaire dans la feuille
  prob<-c()#Proportion de Y=1 dans la feuille
  pop<-list()#Indices des individus dans le noeuds
  tree_split<-list()#Information sur toutes les divisions possibles du noeud
  tables_coupures<-list()#un data.frame pour chaque coupure qui permet de detailler comment l'arbre de coupure est fait
  improvment<-c()#Gain d'impureté dû à la coupure
  action<-c()# 1: si on coupe; -1: quand au moins un critère d'arret rempli; -2: pas d'amelioration de l'impurete qqsoit le groupe ; -3
  parent[1]<-NA
  depth[1]<-0
  groups_selec<-NULL
  n[1]<-dim(data)[1]
  pop[[1]]<-rownames(data)
  i<-1 #indice pour le parcours des feuilles
  i_node_coupure<-NA#indice du noeud dans l'arbre cart de coupure
  i_node_coupure_2<-NA#indice du noeud dans l'arbre cart de coupure
  yval[1] <- ifelse((length(which(data$Y == "1")) / n[1]) < 0.5, "0", "1")
  while(i<=length(n)){
    node<-data[intersect(unlist(pop[[i]]),rownames(data)),]
    prob[i]<-round((length(which(node$Y=="1"))/n[i]),4)
    n_case[i]<-length(which(node$Y=="1"))
    n_noncase[i]<-length(which(node$Y=="0"))
    if(crit=="1"){
      node_impurity<-gini(prop.table(table(node$Y)))
    }else{
      if(crit=="2"){
        node_impurity<-entropy(prop.table(table(node$Y)))
      }else{
        if(crit=="3"){
          node_impurity<-(length(which(label!=node$Y)))- length(which(pred[,j]!=node$Y))
        }
      }
    }
    if ((n_case[i] > case_min) & (n_noncase[i] > case_min)) {### Critères d'arrêt (Nbr min d'observations ou feuille "pure")
      igroups<-group.selection(group[-1],mtry=p)
      groups_selec<-rbind(groups_selec,igroups)
      tree_splits<-list()
      improvment_splits<-rep(NA,length(igroups))
      for(l in 1:length(igroups)){
        tree_splits[[l]]<-cartgv_split(data=node[,c(1,which(group==igroups[l]))],
                                       group=c(NA,1:length(which(group==igroups[l]))),
                                      crit=crit,case_min=1,maxdepth=maxdepth,
                                      p=min(mtry_var[igroups[l]], 
                                            length(which(group==igroups[l]))),penalty="No")
        improvment_splits[l]<-length(node$Y)*(node_impurity-impurity.cartgv(node[,c(1,which(group==igroups[l]))], 
                                                                            list(tree_splits[[l]]$tree),
                                                                            tree_splits[[l]])$impurete[,crit])                       
        group.size<-length(which(group==igroups[l]))
        if(penalty=="Size"){
          improvment_splits[l]<-improvment_splits[l]/group.size
        }
        if(penalty=="Root.size"){
          improvment_splits[l]<-improvment_splits[l]/sqrt(group.size)
        }
        if(penalty=="Log"){
          if(group.size>1){
            improvment_splits[l]<-improvment_splits[l]/log(group.size)
          }
        }
      }
      improvment[i]<-max(improvment_splits)
      if (improvment[i] > 0){#y-a-t-il un groupe de variable qui ameliore l'impurete du noeud?
        action[i] <- 1
        max.importance<-max(improvment_splits)
        if(length(which(improvment_splits==max.importance))>1){
          ind_var<-sample(which(improvment_splits==max.importance),1,FALSE)
        }else{
          ind_var<-which.max(improvment_splits)
        }
        var[i] <- igroups[ind_var]#si oui, on choisit celui qui maximise la decroissance d'impurete
        tree_split[[i]]<-tree_splits[[ind_var]]
        leaves<-which(tree_split[[i]]$tree$leave=="*")
        for(k in 1:length(leaves)){
          i_node_coupure_2<-c(i_node_coupure_2,as.numeric(as.character(tree_split[[i]]$tree$node[leaves[k]])))
          i_node_coupure<-c(i_node_coupure,as.numeric(as.character(tree_split[[i]]$tree$i_node_coupure[leaves[k]])))
          parent<-c(parent,i)
          depth<-c(depth,depth[i]+1)
          yval<-c(yval,as.numeric(as.character(tree_split[[i]]$tree$yval[leaves[k]])))
          pop[[length(pop)+1]]<-tree_split[[i]]$pop[[leaves[k]]]
          n<-c(n,length(unlist(tree_split[[i]]$pop[[leaves[k]]])))
        }
        i<-i+1
        
      }else{
        action[i]<--2 #Aucun groupe de variable n'améliore l'impureté
        var[i]<-NA
        i<-i+1
      }
    }else{
      action[i]<--1 #Au moins un critère d'arrêt rempli
      improvment[i]<-NA
      var[i]<-NA
      i<-i+1
    }
    pred<-NULL
    tree_splits<-NULL
    ind_var<-NULL
    leaves<-NULL
  }
  node<-1:length(depth)
  leave<-ifelse(action<0,"*","")
  improvment<-round(improvment,4)
  tree<-as.data.frame(cbind(action,var,depth,parent,n,n_case,n_noncase,yval,prob,leave,node,improvment,i_node_coupure,i_node_coupure_2))
  rownames(groups_selec)<-paste("split",1:nrow(groups_selec),seq=" ")
  colnames(groups_selec)<-paste(1:ncol(groups_selec),"th group selected",seq=" ")

  return(list(tree=tree,tree_split=tree_split,pop=pop,groups_selec=groups_selec))
}


#' predict.cartgv
#'
#' la fonction fait appel à la fonction predict.cartgv()
#' IMPORTANT : i_noeuds (indice du noeud dans l'arbre de coupure, 
#' càd indice donné par rpart) et noeuds (indice du noeud dans l'arbre cartgv) 
#' utilisées comme clés primaires pour identifier de manière unique un noeuds
#'
#' @param new 
#' @param tree 
#' @param tree_split 
#'
#' @return
#' @export
#'
#' @examples
predict.cartgv.rf<-function(new,tree,tree_split){
  f <- function(x) as.numeric(as.character(x))
  indx <- colnames(tree)
  indx <- indx[which(indx != 'leave')]
  tree[indx] <- lapply(tree[indx], f)
  indx <- colnames(new)
  new[indx]  <- lapply(new[indx],  f)
  
  
  P<-dim(new)[1]
  pred<-rep(NA,P)
  noeuds <- rep(NA, P)
  i_noeuds <- rep(NA, P)
  score <- rep(NA, P)
  for(p in 1:P){
    ind<-new[p,]
    i<-1
    while(tree$action[which(tree$node==i)]>=1){
      temp<-predict.cartgv(ind,
                           tree_split[[i]]$tree,
                           tree_split[[i]]$carts,
                           tree_split[[i]]$tables_coupures)
      
      i<-tree$node[which(tree$parent==i
                           & tree$i_node_coupure_2 == temp$noeuds
                           & tree$i_node_coupure   == temp$i_noeuds)]
                    
    }
    if(tree$action[which(tree$node==i)]<1){
      #on teste si le noeuds dans lequel tombre l'obs est une feuille 
      pred[p]<-tree$yval[which(tree$node==i)]
      noeuds[p] <- i
      i_noeuds[p]<tree$i_node_coupure[which(tree$node==i)]
      score[p]<- tree$prob[which(tree$node==i)]
      i<-1
    }
  }
  
  res<-as.data.frame(cbind(Y=new$Y,
                           hat.Y=pred,
                           noeuds=as.character(noeuds),
                           score=score,
                           i_noeuds=i_noeuds))
  
  names(res)<-c("Y","hat.Y","noeuds","score","i_noeuds")
  if(2%in%res$Y){res$Y<-ifelse(res$Y==2,1,0)}
  if(2%in%res$hat.Y){res$hat.Y<-ifelse(res$hat.Y==2,1,0)}
  return(res)
}





#' Elagage 
#' 
#' calcul de l'impurete 
#' 
#' a partir d'un echantillon independant et d'une sequence d'arbres emboites
#'
#' @param validation 
#' @param tree_seq 
#' @param tree 
#'
#' @return
#' @export
#'
#' @examples
impurity.cartgv.rf <- function(validation, tree_seq,tree) {
  N_tree <- length(tree_seq)
  N <- dim(validation)[1]
  pred <- list()
  sum_noeuds <- list()
  impurete <- matrix(rep(NA, N_tree * 3), nrow = N_tree, ncol = 3)
  for (k in 1:N_tree) {
    predictions <-as.data.frame(predict.cartgv.rf(validation,as.data.frame(tree_seq[[k]]),tree$tree_split))
    predictions<-predictions[order(as.numeric(as.character(predictions$noeuds))),]
    n <- as.numeric(table(as.numeric(as.character(predictions$noeuds))))
    n1 <- tapply(as.numeric(as.character(predictions$Y)), as.numeric(as.character(predictions$noeuds)),sum)
    p1 <- n1 / n
    p0 <- 1 - p1
    predictions$error.pred <-as.numeric(as.character(apply(predictions[, c("hat.Y", "Y")], 1, function(x)ifelse(x[1] != x[2], "1", "0"))))
    misclass <-tapply(as.numeric(as.character(predictions$error.pred)), as.numeric(as.character(predictions$noeuds)), sum)
    summaryNoeuds <-as.data.frame(cbind(as.numeric(as.character(names(table(as.numeric(as.character(predictions$noeuds)))))), n, n1, p1, p0, misclass))
    names(summaryNoeuds) <-c("nom_noeuds","N","N[Y=1]","P[Y=1]","P[Y=0]","P[hat.Y!=Y]")
    impurete[k, 1] <- sum(apply(summaryNoeuds, 1, function(x)x[2] * gini(x[c(4, 5)]))) / N
    impurete[k, 2] <-sum(apply(summaryNoeuds, 1, function(x)x[2] * entropy(x[c(4, 5)]))) / N
    impurete[k, 3] <-sum(summaryNoeuds[, 6]) / (N)
    impurete<-as.data.frame(impurete)
    names(impurete)<-c("Gini","Information","Misclass")
    pred[[k]]<-predictions
    sum_noeuds<-summaryNoeuds
  }
  list(impurete = impurete,pred = pred,summary_noeuds = sum_noeuds)
}

#' grpimpperm.rf
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
grpimpperm.rf<-function(num_group,
                        data,
                        oobsamples,
                        group,
                        tree,
                        impurityacc)
  {
  data_perm<-perm(oobsamples,data,num_group,group)
  accperm<-1-impurity.cartgv.rf(data_perm,list(tree$tree),tree)$impurete$Misclass# accurancy=1-misclass
  DecreaseAcc<-impurityacc-accperm
  return(DecreaseAcc)
}


#' cartgv_split()
#'
#' Utilisé quand sampvar=TRUE et sampvar_type="2" et maxdepth>1 
#' ==> dans l'arbre de coupure avant chaque coupure on tire un sous ensemble de variable
#' Cette fonction faire un arbre de coupure pour un group donné 
#' (donc group=1:nombre de variables dans le group).
#'
#' @param data 
#' @param group 
#' @param crit 
#' @param case_min 
#' @param maxdepth 
#' @param p 
#' @param penalty 
#'
#' @return
#' @export
#'
#' @examples
cartgv_split<-function(data,group,crit=1,case_min=1,maxdepth=2,p=floor(sqrt(length(unique(group[!is.na(group)])))),penalty="No"){
  
  ##Initialisation
  tree<-NULL
  parent<-c()#Noeud ancetre du noeud
  depth<-c()#Profondeur du noeuds
  var<-c()#Groupe selectionné pour couper
  carts<-list()#CART définissant la coupure
  n<-c()# Effectif dans le noeud considéré
  n_case<-c()# Effectif de Y=1 dans le noeud considéré
  n_noncase<-c()# Effectif de Y=0 dans le noeud considéré
  yval<-c()#Prédiction de CART = règle de prediction basée sur la classe majoritaire dans la feuille
  prob<-c()#Proportion de Y=1 dans la feuille
  pop<-list()#Indices des individus dans le noeuds
  splits<-list()#Information sur toutes les divisions possibles du noeud
  tables_coupures<-list()#un data.frame pour chaque coupure qui permet de detailler comment l'arbre de coupure est fait
  improvment<-c()#Gain d'impureté dû à la coupure
  action<-c()# 1: si on coupe; -1: quand au moins un critère d'arret rempli; -2: pas d'amelioration de l'impurete qqsoit le groupe ; -3
  parent[1]<-NA
  depth[1]<-0
  groups_selec<-NULL
  n[1]<-dim(data)[1]
  pop[[1]]<-rownames(data)
  i<-1 #indice pour le parcours des feuilles
  i_node_coupure<-NA#indice du noeud dans l'arbre cart de coupure
  yval[1] <- ifelse((length(which(data$Y == "1")) / n[1]) < 0.5, "0", "1")
  while(i<=length(n)){
    node<-data[intersect(unlist(pop[[i]]),rownames(data)),]
    prob[i]<-round((length(which(node$Y=="1"))/n[i]),4)
    n_case[i]<-length(which(node$Y=="1"))
    n_noncase[i]<-length(which(node$Y=="0"))
    if (((n_case[i] > case_min) & (n_noncase[i] > case_min)) & depth[i]<maxdepth) {### Critères d'arrêt (Nbr min d'observations ou feuille "pure")
      igroups<-group.selection(group[!is.na(group)],mtry=p)
      groups_selec<-rbind(groups_selec,igroups)
      splits[[i]]<-split.cartgv(node=node,group=group,igroups,label=yval[i],maxdepth=1,penalty="No",sampvar=FALSE)
      improvment[i]<-max(unlist(splits[[i]][[crit]]))
      if (improvment[i] > 0){#y-a-t-il un groupe de variable qui ameliore l'impurete du noeud?
        action[i] <- 1
        max.importance<-max(unlist(splits[[i]][[crit]]))
        if(length(which(unlist(splits[[i]][[crit]])==max.importance))>1){
          ind_var<-sample(which(unlist(splits[[i]][[crit]])==max.importance),1,FALSE)
        }else{
          ind_var<-which.max(unlist(splits[[i]][[crit]]))
        }
        var[i] <- igroups[ind_var]#si oui, on choisit celui qui maximise la decroissance d'impurete
        carts[[i]]<-splits[[i]]$carts[[ind_var]]
        tables_coupures[[i]]<-calcul_cart(carts[[i]],node)
        modalities<-unique(as.numeric(as.character(carts[[i]]$where)))
        for(k in 1:length(modalities)){
          i_node_coupure<-c(i_node_coupure,modalities[k])
          node_son<-node[which(unlist(carts[[i]]$where)==modalities[k]),]
          parent<-c(parent,i)
          depth<-c(depth,depth[i]+1)
          n<-c(n,dim(node_son)[1])
          yval<-c(yval,unique(splits[[i]]$pred[which(unlist(carts[[i]]$where)==modalities[k]),ind_var]))
          pop[[length(pop)+1]]<-rownames(node_son)
        }
        i<-i+1
        
      }else{
        action[i]<--2 #Aucun groupe de variable n'améliore l'impureté
        var[i]<-NA
        i<-i+1
      }
    }else{
      action[i]<--1 #Au moins un critère d'arrêt rempli
      improvment[i]<-NA
      var[i]<-NA
      i<-i+1
    }
    pred<-NULL
    cart<-NULL
    ind_var<-NULL
    modalities<-NULL
  }
  node<-1:length(depth)
  leave<-ifelse(action<0,"*",ifelse(depth==maxdepth,"*",""))
  improvment<-round(improvment,4)
  tree<-as.data.frame(cbind(action,var,depth,parent,n,n_case,
                            n_noncase,yval,prob,leave,node,
                            improvment,i_node_coupure))
  rownames(groups_selec)<-paste("split",1:nrow(groups_selec),seq=" ")
  colnames(groups_selec)<-paste(1:ncol(groups_selec),"th group selected",seq=" ")
  return(list(tree=tree,carts=carts,splits=splits,pop=pop,
              tables_coupures=tables_coupures,groups_selec=groups_selec))
}



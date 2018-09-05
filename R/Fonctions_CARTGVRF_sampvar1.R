
as.numeric.factor <- function(x) {as.numeric(as.character(x))}

#' cartgv
#' 
#' Classification And Regression Trees for Grouped Variables
#' 
#' Function which takes as inputs a integer "B"  and a sample of data "data" . The function 
#' draws randomly B bootstrap samples of size N (the size of the sample of data) and returns 
#' a matrix which with N lines and B columns which contains the indices of the observations 
#' belonging each boostrap sample. 
#' 
#' @param data the trainning dataset
#' @param group  a vector with the group number of each variable (WARNING : if there are P goups, then the vector group must contains all the numbers from 1 to P )
#' @param crit 1=Gini,2=Entropie,3=Misclassification error
#' @param case_min minimum number of cases/non cases in a node which allows to split the node
#' @param maxdepth the max depth for a split-tree
#' @param RF 
#' @param penalty "No","Size","Root.size" or "Log"
#' @param p 
#' @param IMPORTANCE 
#' @param sampvar a boolean indicating if within each tree-split, a subset of variables is drawn for each group
#' @param mtry_var usefull only if sampvar=TRUE. It indicates the number of drawn variables
#'
#' @return
#'    - tree : data.frame "resumant" l'arbre construit 
#'    - carts : liste des arbres carts utilises pour elaborer l'arbre
#'    - splits : liste donnant pour chaque noeud la valeur de la fonction Split_Tree_CART()
#'    - pop : liste contenant pour chaque neoud le nom des individus (rownames)
#'    - tables_coupures : liste contenant pour chaque coupure un data.frame permettant de reconstruire l'arbre CART de coupure
#'    - groups_selec   : liste contenant pour chaque coupure les groupes candidats
#'    - importance     : liste donnant pour chaque groupe la décroissance d'impureté resultant de la meilleur coupure opérée sur le groupe
#'    - agreement1      : liste donnant pour chaque coupure et chaque groupe candidat la proba d'accord (Indice de Rand) entre la coupure choisie et la meilleure coupure pour le groupe
#'    - cimportance1     : liste donnant pour chaque coupure et chaque groupe l'importance "corrigée" par la proba d'accord
#' @export
#'
#' @examples
cartgv<-function(data,
                 group,
                 crit=1,
                 case_min=1,
                 maxdepth=2,
                 RF=FALSE,
                 penalty="No",
                 p=ifelse(!is.null(mtry_group),mtry_group,sqrt(length(unique(group[!is.na(group)])))), 
                 IMPORTANCE=ifelse(RF==TRUE,FALSE,TRUE),
                 sampvar=ifelse(RF==TRUE,TRUE,FALSE),
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
  splits<-list()#Information sur toutes les divisions possibles du noeud
  tables_coupures<-list()#un data.frame pour chaque coupure qui permet de detailler comment l'arbre de coupure est fait
  improvment<-c()#Gain d'impureté dû à la coupure
  action<-c()# 1: si on coupe; -1: quand au moins un critère d'arret rempli; -2: pas d'amelioration de l'impurete qqsoit le groupe ; -3
  if(IMPORTANCE==T){
    importance<-list()
    cimportance<-list()
    agreement<-list()
  }
  parent[1]<-NA
  depth[1]<-0
  groups_selec<-NULL
  n[1]<-dim(data)[1]
  pop[[1]]<-rownames(data)
  i<-1 #indice pour le parcours des feuilles
  i_node_coupure<-NA#indice du noeud dans l'arbre cart de coupure
  yval[1] <- ifelse((length(which(data$Y == 1)) / n[1]) < 0.5, 0, 1)
  while(i<=length(n)){
    node<-data[intersect(unlist(pop[[i]]),rownames(data)),]
    prob[i]<-round((length(which(node$Y==1))/n[i]),4)
    n_case[i]<-length(which(node$Y==1))
    n_noncase[i]<-length(which(node$Y==0))
    if ((n_case[i] > case_min) & (n_noncase[i] > case_min)) {### Critères d'arrêt (Nbr min d'observations ou feuille "pure")
      if(RF==TRUE){
        igroups<-group.selection(group[!is.na(group)],mtry=p)
        groups_selec<-rbind(groups_selec,igroups)
      }else{
        igroups<-unique(group[!is.na(group)])
      }
      splits[[i]]<-split.cartgv(node=node,group=group,igroups,label=yval[i],maxdepth=maxdepth,penalty=penalty,sampvar=sampvar,mtry_var=mtry_var)
      improvment[i]<-max(unlist(splits[[i]][[crit]]))
      if (improvment[i] > 0){#y-a-t-il un groupe de variable qui ameliore l'impurete du noeud?
        action[i] <- 1
        max.importance<-max(unlist(splits[[i]][[crit]]))
        if(length(which(unlist(splits[[i]][[crit]])==max.importance))>1){
          ind_var<-sample(which(unlist(splits[[i]][[crit]])==max.importance),1,FALSE)
        }else{
          ind_var<-which.max(unlist(splits[[i]][[crit]]))
        }
        
        
        if(IMPORTANCE==TRUE){
          importance[[i]]<-unlist(splits[[i]][[crit]])
          agreement[[i]]<-sapply(1:length(igroups),
                                  function(x) e1071::classAgreement(table(unlist(splits[[i]]$carts[[x]]$where),
                                                                         unlist(splits[[i]]$carts[[ind_var]]$where)),match.names=FALSE)$rand)
          cimportance[[i]]<-unlist(importance[[i]])*unlist(agreement[[i]])
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
  leave<-ifelse(action<0,"*","")
  improvment<-round(improvment,4)
  tree<-as.data.frame(cbind(action,
                            var,
                            depth,
                            parent,
                            n,n_case,
                            n_noncase,
                            yval,
                            prob,
                            leave,
                            node,
                            improvment,
                            i_node_coupure))
  if(RF==T){
    rownames(groups_selec)<-paste("split",1:nrow(groups_selec),seq=" ")
    colnames(groups_selec)<-paste(1:ncol(groups_selec),"th group selected",seq=" ")
  }
  if(IMPORTANCE==FALSE){
    importance<-NULL
    cimportance<-NULL
    agreement<-NULL
  }
  return(list(tree=tree,
              carts=carts,
              splits=splits,
              pop=pop,
              tables_coupures=tables_coupures,
              groups_selec=groups_selec,
              importance=importance,
              cimportance=cimportance,
              agreement=agreement))
}

#' split.cartgv()
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
split.cartgv<-function(node,
                       group,
                       igroups,
                       label,
                       maxdepth=2,
                       penalty="No",
                       sampvar="FALSE",
                       mtry_var=sapply(as.numeric(table(group[!is.na(group)])),function(x)floor(sqrt(x))))
  {
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
      ivar<-sample(which(group==igroups[j]),
                   size=min(mtry_var[igroups[j]], 
                            length(which(group==igroups[j]))),replace=FALSE)
      var_selec[ivar]<-1
    }else{
      ivar<-which(group==igroups[j])
    }
    formula<-stats::as.formula(paste("Y ~ ",paste(names(node)[ivar],collapse="+")))
    cart <- rpart(formula,
                  node,
                  control=rpart.control(minsplit=1,
                                        cp=-0.999,
                                        maxdepth = maxdepth,
                                        maxsurrogate =0,
                                        maxcompete = 0 ),
                  method = "class",
                  parms = list(split = "gini"))
    
    pred[,j]<-as.numeric(as.character(stats::predict(cart, type = "class", newdata=node)))
    carts[[j]] <- cart
    cart<-NULL
  }
  
  for(j in 1:nb_group){
    if(sum(pred[,j]==2,na.rm=T)>0){pred[,j]<-ifelse(pred[,j]==2,1,0)}
    tab<-table(pred[,j],node$Y)
    modalities<-unique(carts[[j]]$where)
    Gain_Gini[j]<-length(node$Y)*(gini(prop.table(table(node$Y))))
    Gain_Ent[j]<-length(node$Y)*(entropy(prop.table(table(node$Y))))
    Gain_Clas[j]<-length(node$Y)*((length(which(label!=node$Y)))- length(which(pred[,j]!=node$Y)))
    
    for(k in 1:length(modalities)){
      
      node1<-node[which(unlist(carts[[j]]$where)==modalities[k]),]
      tab1<-table(pred[which(unlist(carts[[j]]$where)==modalities[k]),j],
                  node$Y[which(unlist(carts[[j]]$where)==modalities[k])])
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
  return(list(Gain_Gini=Gain_Gini, 
              Gain_Ent=Gain_Ent, 
              Gain_Clas=Gain_Clas, 
              carts=carts, 
              pred=pred,
              igroups=igroups,
              var_selec=var_selec))
}


#' perm
#'
#' @param oobsamples 
#' @param data 
#' @param num_group 
#' @param group 
#'
#' @return
#'
#' @examples
perm<-function(oobsamples,data,num_group,group){
  data_perm<-data[oobsamples,]
  iperm<-gtools::permute(1:nrow(data_perm))
  data_perm[,which(group==num_group)]<-data_perm[iperm,which(group==num_group)]
  names(data_perm)<-names(data)
  return(data_perm)
}


#' grpimpperm
#'
#' @param num_group 
#' @param data 
#' @param oobsamples 
#' @param group 
#' @param tree 
#' @param impurityacc 
#'
#' @examples
grpimpperm<-function(num_group,data,oobsamples,group,tree,impurityacc){
  data_perm<-perm(oobsamples,data,num_group,group)
  accperm<-1-impurity.cartgv(data_perm,list(tree$tree),tree)$impurete$Misclass# accurancy=1-misclass
  DecreaseAcc<-impurityacc-accperm
  return(DecreaseAcc)
}


#' grpimpgini
#'
#' @param num_group 
#' @param groupselec 
#' @param tree 
#'
#' @return
#'
#' @examples
grpimpgini<-function(num_group,groupselec,tree){
  
  
  DecreaseImpurity<-0
  times<-which(tree$tree$var==num_group)
  if(num_group%in%groupselec){
    for(i in times ){
      p<-unlist(as.numeric.factor(tree$tree$n[i]))/unlist(as.numeric.factor(tree$tree$n[1]))
      DecreaseImpurity<-DecreaseImpurity+(p*unlist(tree$split[[i]]$Gain_Gini[which(tree$split[[i]]$igroups==num_group)]))
    }
  }
  return(DecreaseImpurity)
}


#' predict.cartgv
#' 
#' Predict classification and regression tree for grouped variables
#' 
#'la fonction fait appel a la fonction pred_cart()
#' IMPORTANT : i_noeuds (indice du noeud dans l'arbre de coupure, cad indice donne par rpart) 
#' et noeuds (indice du noeud dans l'arbre cartgv) utilisees comme cles primaires pour 
#' identifier de manière unique un noeuds
#'
#' @param new 
#' @param tree 
#' @param carts 
#' @param coups 
#'
#' @return
#'  Valeurs retournees par la fonction : 
#' la fonction renvoie une matrice donnant pour chaque individu de new (dans l'ordre) : 
#'  hat.Y :la classe predite
#'  Y : la classe d'appartenance
#'  noeuds : le numero du noeud contenant l'observation
#'  score : le score le l'observation

#' @export
#'
#' @examples
predict.cartgv<-function(new,tree,carts,coups){
  
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
      temp<-as.data.frame(t(pred_cart(ind,coups[[i]])[c(2,3,4,5)]))
      names(temp)<-c("i_node_coupure","pred","n_noncase","n_case")
      i<-tree$node[which(tree$parent==i
                        & tree$yval==temp$pred
                        & tree$n_case==temp$n_case
                        & tree$n_noncase==temp$n_noncase
                        & tree$i_node_coupure== temp$i_node_coupure)]
    }
    if(tree$action[which(tree$node==i)]<1){
      #on teste si le noeuds dans lequel tombre l'obs est une feuille 
      pred[p]<-tree$yval[which(tree$node==i)]
      noeuds[p] <- i
      i_noeuds[p]<-tree$i_node_coupure[which(tree$node==i)]
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


#' Title
#'
#' @param new 
#' @param tree 
#' @param carts 
#' @param coups 
#'
#' @return
#' la fonction renvoie une matrice donnant pour chaque 
#' individu de new (dans l'ordre) :
#'  - hat.Y :la classe predite
#'  - Y : la classe d'appartenance
#'  - noeuds : le numero du noeud contenant l'observation
#'  - score : le score le l'observation
#'
#' @examples
predict.test.cartgv <- function(new, tree, carts, coups) {
  P <- dim(new)[1]
  
  predgv <- function(ind) {
    i <- 1
    while (tree$action[which(tree$node == i)] >= 1) {
      tab <- coups[[i]]
      if (is.data.frame(tab)) {
        i_node <- 1
        while (tab$var[i_node] != "<leaf>") {
          i_sons <- which(tab$parent == tab$node_name_cart[i_node])
          if (ind[tab$var[i_node]] >= tab$threshold[i_node]) {
            if (tab$sens[i_node] == "sup")
              k <- 1
            else
              k <- 2
          } else{
            if (tab$sens[i_node] == "sup")
              k <- 2
            else
              k <- 1
          }
          i_node <- i_sons[k]
        }
        
        i <-  tree$node[which(
          tree$parent == i
          & tree$i_node_coupure == i_node
          & tree$yval == tab$pred[i_node]
          & tree$n_noncase == tab$n_noncase[i_node]
          & tree$n_case == tab$n_case[i_node]
        )]
      }
    }
    i.eq.tree.node <- which(tree$node == i)
    if (tree$action[i.eq.tree.node] < 1) {
      # on teste si le noeuds dans lequel tombre l'obs
      # est une feuille
      return(c(tree$yval[i.eq.tree.node], i, tree$prob[i.eq.tree.node]))
    } else {
      return(c(NA, NA, NA))
    }
    
  }
  
  res <- matrix()
  length(res) <- 4 * P
  dim(res) <- c(P, 4)
  
  res[,1] <- new[,1]
  res[,2:4] <- t(apply(new,1,predgv))
  
  colnames(res) <- c("Y", "hat.Y", "noeuds", "score")
  
  return(as.data.frame(res))
}



#' Title
#' Elagage :  calcul de l'impuret'e a' partir d'un 'echantillon independant et d'une s'equence d'arbres emboites
#' WARNING : fonction corrigée le 19 janvier2018 car il y avait un problème d'ordre pour les noeuds dans les apply
#' La fonction prend en entrée un échantillon test et une séquence d'arbres emboités
#' Pour chaque arbre, on prédit la classe de chaque observation de l'ensemble test;
#' puis à partir de ces résultats, on calcul l'impureté de chaque arbre (gini, l'entropie et le taux de mal-classés sont calculés)
#' method : la nature des coupures pour l'arbre CC binaire ("bin" si tree construit avec Tree_CART_bin()) 
#'          ou ("multi" si tree construit avec Tree_CART()) 
#'
#' @param validation data.frame contenant les mêmes variables que le jeux de données utilisés pour construire l'objet tree
#' @param tree_seq objet retourné par la fonction Tree_CART_bin() ou Tree_CART()
#' @param tree 
#'
#' @return
#' List : 
#'   impurete : une matrice contenant les differentes valeurs d'impuret'e pour chaque sous-arbres
#'   pred = liste des predictions pour chaque sous arbres
#'   summary_noeuds : liste contenant pour chaque sous-arbre des infos sur ses noeuds: 
#'                     nom_noeuds = noms du noeuds
#'                     N = nb obs dans le neoud
#'                     N[Y=1] = nb d'obs dans le noeuds appartenant a' la class "Y=1" 
#'                     P[Y=1] = proba pour une obs du noeuds d'appartenanir  a' la class "Y=1"
#'                     P[Y=0] = proba pour une obs du noeuds d'appartenanir  a' la class "Y=0"
#'                     P[hat.Y!=Y] = tx de mal class'es dans le noeud
#' @export
#'
#' @examples
impurity.cartgv <- function(validation, tree_seq,tree) {
  N_tree <- length(tree_seq)
  N <- dim(validation)[1]
  pred <- list()
  sum_noeuds <- list()
  impurete <- matrix(rep(NA, N_tree * 3), nrow = N_tree, ncol = 3)
  for (k in 1:N_tree) {
    predictions <-as.data.frame(predict.cartgv(validation,
                                               as.data.frame(tree_seq[[k]]),
                                               tree$carts,
                                               tree$tables_coupures))
    
    predictions<-predictions[order(as.numeric(as.character(predictions$noeuds))),]
    n <- as.numeric(table(as.numeric(as.character(predictions$noeuds))))
    n1 <- tapply(as.numeric.factor(predictions$Y), as.numeric.factor(predictions$noeuds),sum)
    p1 <- n1 / n
    p0 <- 1 - p1
    predictions$error.pred <-as.numeric.factor(apply(predictions[, c("hat.Y", "Y")], 
                                                           1, 
                                                           function(x)ifelse(x[1] != x[2], 1, 0)))
    misclass <-tapply(as.numeric.factor(predictions$error.pred), 
                      as.numeric.factor(predictions$noeuds), 
                      sum)
    summaryNoeuds <-as.data.frame(cbind(as.numeric.factor(names(table(as.numeric.factor(predictions$noeuds)))), 
                                        n, n1, p1, p0, misclass))
    names(summaryNoeuds) <-c("nom_noeuds","N","N[Y=1]","P[Y=1]","P[Y=0]","hat.Y!=Y")
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


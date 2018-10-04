############################################################################################
#################           Functions about the variant of              ####################
#################        the CARTGV approach for RFGV forests           ####################
############################################################################################
as.numeric.factor <- function(x) {as.numeric(as.character(x))}

#' cartgv.rf
#'
#' Variant of the CARTGV approach to build RFGV forests. Implemented for binary classification problems
#'
#' @param data a data frame containing the response value (for the first variable)  and the predictors and used to grow the tree.
#' The name of the response value must be "Y".The response variable must be the first variable of the data frame and the variable
#' must be coded as the two levels "0" and "1".
#' @param group  group a vector with the group number of each variable.
#' (WARNING : if there are "\code{p}" goups, the groups must be numbers from "\code{1}" to "\code{p}" in increasing order. The
#' group label of the response variable is missing (i.e. NA)).
#' @param crit an integer indicating the impurity function used (1=Gini index / 2=Entropie/ 3=Misclassification rate).
#' @param penalty a boolean indicating if the decrease in node impurity must take account of the group size. Four penalty are
#' available: "No","Size","Root.size" or "Log".
#' @param case_min an integer indicating the minimun number of cases/non cases in a terminal nodes. The default is 1.
#' @param maxdepth the max depth for a split-tree.
#' @param mtry_group an integer the number of variables randomly samples as candidates at each split.
#' @param mtry_var a vector of length the number of groups. It indicates the number of drawn variables for each group.

#' @return a list with elements
#'    - tree: a data frame which summarizes the resulted CARTGV tree.
#'    - tree_split: a list containing informations about the splitting trees. Each element is an object returned by the function "\code{cartgv_split}".
#'    - pop: a list containing the indices (rownames) of the observations which belong to the nodes.
#'    - groups_selec: a matrix containint for each splitting-tree the indices of the sampled grouped. Precisely, the i-th row correspond to the i-th splitting-tree.
#
#' @export
#'
cartgv.rf<-function(data,
                    group,
                    crit=1,
                    case_min=1,
                    maxdepth=2,
                    mtry_group=floor(sqrt(length(unique(group[!is.na(group)])))),
                    penalty="No",
                    mtry_var=sapply(as.numeric(table(group[!is.na(group)])),function(x)floor(sqrt(x))))
  {

  tree<-NULL
  parent<-c()
  depth<-c()
  var<-c()
  carts<-list()
  n<-c()
  n_case<-c()
  n_noncase<-c()
  yval<-c()
  prob<-c()
  pop<-list()
  tree_split<-list()
  tables_coupures<-list()
  improvment<-c()
  action<-c()
  parent[1]<-NA
  depth[1]<-0
  groups_selec<-NULL
  n[1]<-dim(data)[1]
  pop[[1]]<-rownames(data)
  i<-1
  i_node_coupure<-NA#node index in the CART tree used to split the node into two new nodes
  i_node_coupure_2<-NA#node index in the splitting tree
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
          node_impurity<-length(which(yval[i]!=node$Y))
        }
      }
    }
    if ((n_case[i] > case_min) & (n_noncase[i] > case_min)) {
      igroups<-group.selection(group[-1],mtry=mtry_group)
      groups_selec<-rbind(groups_selec,igroups)
      tree_splits<-list()
      improvment_splits<-rep(NA,length(igroups))

      for(l in 1:length(igroups)){
        tree_splits[[l]]<-cartgv_split(data=node[,c(1,which(group==igroups[l]))],
                                       group=c(NA,1:length(which(group==igroups[l]))),
                                       crit=crit,
                                       case_min=1,
                                       maxdepth=maxdepth,
                                       p=min(mtry_var[igroups[l]],
                                            length(which(group==igroups[l]))),
                                       penalty="No")

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
      if (improvment[i] > 0){
        action[i] <- 1
        max.importance<-max(improvment_splits)
        if(length(which(improvment_splits==max.importance))>1){
          ind_var<-sample(which(improvment_splits==max.importance),1,FALSE)
        }else{
          ind_var<-which.max(improvment_splits)
        }
        var[i] <- igroups[ind_var]
        tree_split[[i]]<-tree_splits[[ind_var]]
        leaves<-which(tree_split[[i]]$tree$leave=="*")
        for(k in 1:length(leaves)){
          i_node_coupure_2<-c(i_node_coupure_2,as.numeric.factor(tree_split[[i]]$tree$node[leaves[k]]))
          i_node_coupure<-c(i_node_coupure,as.numeric.factor(tree_split[[i]]$tree$i_node_coupure[leaves[k]]))
          parent<-c(parent,i)
          depth<-c(depth,depth[i]+1)
          yval<-c(yval,as.numeric.factor(tree_split[[i]]$tree$yval[leaves[k]]))
          pop[[length(pop)+1]]<-tree_split[[i]]$pop[[leaves[k]]]
          n<-c(n,length(unlist(tree_split[[i]]$pop[[leaves[k]]])))
        }
        i<-i+1

      }else{
        action[i]<--2
        var[i]<-NA
        i<-i+1
      }
    }else{
      action[i]<--1
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

  tree<-as.data.frame(cbind(action,
                            var,
                            depth,
                            parent,
                            n,
                            n_case,
                            n_noncase,
                            yval,
                            prob,
                            leave,
                            node,
                            improvment,
                            i_node_coupure,
                            i_node_coupure_2))
  rownames(groups_selec)<-paste("split",1:nrow(groups_selec),seq=" ")
  colnames(groups_selec)<-paste(1:ncol(groups_selec),"th group selected",seq=" ")
  return(list(tree=tree,tree_split=tree_split,pop=pop,groups_selec=groups_selec))
}

#' cartgv_split
#'
#' Build a splitting tree in the modified CARTGV trees containing in a RFGV forest.
#'
#' @param data a data frame containing the response value (for the first variable)
#' and the predictors and used to grow the tree.
#' The name of the response value must be "Y".The response variable must be the
#' first variable of the data frame and the variable
#' must be coded as the two levels "0" and "1".
#' @param group  group a vector with the group number of each variable.
#'  (WARNING : if there are "\code{p}" goups, the groups must be numbers from
#'  "\code{1}" to "\code{p}" in increasing order. The
#' group label of the response variable is missing (i.e. NA)).
#' @param crit an integer indicating the impurity function used
#' (1=Gini index / 2=Entropie/ 3=Misclassification rate).
#' @param penalty a boolean indicating if the decrease in node impurity must take
#' account of the group size. Four penalty are
#' available: "No","Size","Root.size" or "Log".
#' @param case_min an integer indicating the minimun number of cases/non cases in
#' a terminal nodes. The default is 1.
#' @param maxdepth the max depth for a split-tree.
#' @param p an integer indicating the number of variables randomly samples as
#' candidates at each split.
#'
#' @return a list with elements
#'    - tree : a data frame which summarizes the resulted splitting tree.
#'    - carts : a list containing all the CART objects used to buid the splitting
#'    tree. (Note that each split in the splitting tree is a CART object)
#'    - splits : a list containing informations about the splits. Each element
#'    is an object retuned by the function "\code{split_cartgv}".
#'    - pop : a list containing the indices (rownames) of the observations which
#'    belong to the nodes.
#'    - tables_coupures : a list containing data frames that summarizes the splits.
#'    - groups_selec : a matrix containint for each splitting-tree the indices of
#'    the sampled grouped. Precisely, the i-th row correspond to the i-th splitting-tree.
#'
#' @export
cartgv_split<-function(data,
                       group,
                       crit=1,
                       case_min=1,
                       maxdepth=2,
                       p=floor(sqrt(length(unique(group[!is.na(group)])))),
                       penalty="No"){

  ##Initialisation
  tree<-NULL
  parent<-c()
  depth<-c()
  var<-c()
  carts<-list()
  n<-c()
  n_case<-c()
  n_noncase<-c()
  yval<-c()
  prob<-c()
  pop<-list()
  splits<-list()
  tables_coupures<-list()
  improvment<-c()
  action<-c()
  parent[1]<-NA
  depth[1]<-0
  groups_selec<-NULL
  n[1]<-dim(data)[1]
  pop[[1]]<-rownames(data)
  i<-1
  i_node_coupure<-NA
  yval[1] <- ifelse((length(which(data$Y == "1")) / n[1]) < 0.5, "0", "1")
  while(i<=length(n)){
    node<-data[intersect(unlist(pop[[i]]),rownames(data)),]
    prob[i]<-round((length(which(node$Y=="1"))/n[i]),4)
    n_case[i]<-length(which(node$Y=="1"))
    n_noncase[i]<-length(which(node$Y=="0"))
    if (((n_case[i] > case_min) & (n_noncase[i] > case_min)) & depth[i]<maxdepth) {
      igroups<-group.selection(group[!is.na(group)],mtry=p)
      groups_selec<-rbind(groups_selec,igroups)

      splits[[i]]<-split_cartgv(node=node,
                                group=group,
                                label=yval[i],
                                maxdepth=1,
                                penalty="No")

      improvment[i]<-max(unlist(splits[[i]][[crit]]))
      if (improvment[i] > 0){
        action[i] <- 1
        max.importance<-max(unlist(splits[[i]][[crit]]))
        if(length(which(unlist(splits[[i]][[crit]])==max.importance))>1){
          ind_var<-sample(which(unlist(splits[[i]][[crit]])==max.importance),1,FALSE)
        }else{
          ind_var<-which.max(unlist(splits[[i]][[crit]]))
        }
        var[i] <- igroups[ind_var]
        carts[[i]]<-splits[[i]]$carts[[ind_var]]
        tables_coupures[[i]]<-calcul_cart(carts[[i]],node)
        modalities<-unique(as.numeric.factor(carts[[i]]$where))
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
        action[i]<--2
        var[i]<-NA
        i<-i+1
      }
    }else{
      action[i]<--1
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
  return(list(tree=tree,
              carts=carts,
              splits=splits,
              pop=pop,
              tables_coupures=tables_coupures,
              groups_selec=groups_selec))
}


#' predict_cartgv.rf
#'
#' Prediction of test data from a fitted modified CARTGV tree.
#'
#' The function called the function "\code{predict_cartgv}".
#'
#' @param new  a new data frame containing the same variables that "\code{data}".
#' @param tree the data frame "\code{tree}" returned by the function "\code{cartgv.rf}".
#' @param tree_split the object "\code{tree_split}" returned by the function "\code{cartgv.rf}".
#'
#' @return a matrix with "\code{nrows(new)}" rows (the i-th row is provided prediction information about the -ith observation in "\code{new}")
#'         and the 4 followingcolumns,
#'         - Y: the true label ("0" or "1"),
#'         - hat.Y: the predicted label ("0" or "1"),
#'         - noeuds: name of the node in the modified CARTGV tree,
#'         - score:  conditionnal probability that the observation has "Y=1",
#'         - i_noeuds: name of the node in the splitting tree.
#'
#' @export
#'
predict_cartgv.rf<-function(new,tree,tree_split){
  indx <- colnames(tree)
  indx <- indx[which(indx != 'leave')]
  tree[indx] <- lapply(tree[indx], as.numeric.factor)
  indx <- colnames(new)
  new[indx]  <- lapply(new[indx],  as.numeric.factor)

  P<-dim(new)[1]
  pred<-rep(NA,P)
  noeuds <- rep(NA, P)
  i_noeuds <- rep(NA, P)
  score <- rep(NA, P)
  for(p in 1:P){
    ind<-new[p,]
    i<-1
    while(tree$action[which(tree$node==i)]>=1){
      temp<-predict_cartgv(ind,
                           tree_split[[i]]$tree,
                           tree_split[[i]]$carts,
                           tree_split[[i]]$tables_coupures)

      i<-tree$node[which(tree$parent==i
                           & tree$i_node_coupure_2 == temp$noeuds
                           & tree$i_node_coupure   == temp$i_noeuds)]

    }
    if(tree$action[which(tree$node==i)]<1){
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


#' impurity.cartgv.rf
#'
#' Compute the impurity information from a sequence of subtrees by using an independent set.
#' The subtrees are based on a modified CARTGV tree (i.e. a CARTGV tree that is included in a RFGV forest).
#'
#' @param validation a new data frame containing the same variables that "\code{data}".
#' @param tree_seq the object returned by the function "\code{extract_subtrees}". Each element of this sequence is an object "\code{tree}" returned by the function \code{cartgv.rf}.
#' @param tree a fitted modified CARTGV tree. It is an output of the function "\code{cartgv.rf}".
#'
#' @return a list with elements
#'        - impurete: a data frame containing the value of several impurity fucntions (in this order Gini, Entropy, misclassification rate)
#'                    for each subtree of the sequence. The i-th row corresponds to the i-th subtree of the sequence.
#'        - pred:  a list containing the prediction of the label for the data set "\code{validation}" based on each subtree.
#'                 Precisely, the i-th element is the object returned by the function "\code{predict_cartgv.rf}" for the i-th subtree nd by using
#'                 the data set "\code{validation}".
#'        - summary_noeuds: a list containg for each subtree informations about the nodes  (nom_noeuds: node name, N: number of observations in
#'                          the node, \code{N[Y=1]}: number of observation with "Y=1" in the node, \code{N[Y=0]}: number of observation with "Y=0"
#'                          in the node, \code{P[Y=1]}: estimated probability that an observation in the node is assigned to the label "Y=1",
#'                          \code{P[Y=1]}: estimated probability that an observation in the node is assigned to the label "Y=0" and \code{P[hat.Y!=Y]}:
#'                          misclassification rate in the node).
#' @export
#'
impurity.cartgv.rf <- function(validation, tree_seq,tree) {
  N_tree <- length(tree_seq)
  N <- dim(validation)[1]
  pred <- list()
  sum_noeuds <- list()
  impurete <- matrix(rep(NA, N_tree * 3), nrow = N_tree, ncol = 3)
  for (k in 1:N_tree) {
    predictions <-as.data.frame(predict_cartgv.rf(validation,as.data.frame(tree_seq[[k]]),tree$tree_split))
    predictions<-predictions[order(as.numeric.factor(predictions$noeuds)),]
    n <- as.numeric(table(as.numeric.factor(predictions$noeuds)))
    n1 <- tapply(as.numeric.factor(predictions$Y), as.numeric.factor(predictions$noeuds),sum)
    p1 <- n1 / n
    p0 <- 1 - p1
    predictions$error.pred <-as.numeric.factor(apply(predictions[, c("hat.Y", "Y")],
                                                     1, function(x)ifelse(x[1] != x[2], "1", "0")))

    misclass <-tapply(as.numeric.factor(predictions$error.pred),
                      as.numeric.factor(predictions$noeuds), sum)

    summaryNoeuds <-as.data.frame(cbind(as.numeric.factor(names(table(as.numeric.factor(predictions$noeuds)))),
                                        n, n1, p1, p0, misclass))

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
#' Compute the permutation importance of a group from a modified CARTGV tree. Used only by the function \code{rfgv} when grp.importance=TRUE.
#'
#' @param num_group index of the considered group
#' @param data a data frame containing all observation included in the training data set.
#' @param group a vector with the group number of each variable. (WARNING : if there are "\code{p}" goups, the groups must be numbers
#' from "\code{1}" to "\code{p}" in increasing order. The group label of the response variable is missing (i.e. NA)).
#' @param oobsamples a vector containing the indices of the observations that are not included in the sample used to build the modified CARTGV tree.
#' @param tree the output of the function \code{cartgv.rf}.
#' @param impurityacc the accurancy value of the modified CARTGV tree.
#'
#' @return DecreaseAcc the decrease in accurance from permuting the values of the considered group.
#'
grpimpperm.rf<-function(num_group,
                        data,
                        oobsamples,
                        group,
                        tree,
                        impurityacc)
  {
  data_perm<-perm(oobsamples,data,num_group,group)
  accperm<-1-impurity.cartgv.rf(data_perm,
                                list(tree$tree),
                                tree)$impurete$Misclass# accurancy=1-misclass
  DecreaseAcc<-impurityacc-accperm
  return(DecreaseAcc)
}



#' perm
#'
#' Function that permutate the values of a group of variables.
#' Note that all the variables of the group are permuted together in order to keep the relationship between variables within the group.
#'
#' @param oobsamples a vector containing the indices of the observations that are not included in the sample used to build the modified CARTGV tree.
#' @param data a data frame containing all observation included in the training data set.
#' @param num_group index of the considered group
#' @param group a vector with the group number of each variable. (WARNING : if there are "\code{p}" goups, the groups must be numbers
#' from "\code{1}" to "\code{p}" in increasing order. The group label of the response variable is missing (i.e. NA)).
#'
#' @return data_perm a data frame obtained after permuting the values of the considered group for observation belonging to the out-of-bag sample.
#'
perm<-function(oobsamples,data,num_group,group){
  data_perm<-data[oobsamples,]
  iperm<-gtools::permute(1:nrow(data_perm))
  data_perm[,which(group==num_group)]<-data_perm[iperm,which(group==num_group)]
  names(data_perm)<-names(data)
  return(data_perm)
}








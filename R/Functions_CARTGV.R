############################################################################################
#################             Function to build CARTGV trees            ####################
############################################################################################

as.numeric.factor <- function(x) {as.numeric(as.character(x))}

#' cartgv
#'
#' Classification And Regression Trees for Grouped Variables.
#'
#' Implemented for binary classification problems
#'
#' @param data a data frame containing the response value (for the first variable)  and the predictors and used to grow the tree.
#' The name of the response value must be "Y".The response variable must be the first variable of the data frame and the variable
#' must be coded as the two levels "0" and "1".
#' @param group  group a vector with the group number of each variable.
#'  (WARNING : if there are "\code{p}" goups, the groups must be numbers from "\code{1}" to "\code{p}" in increasing order. The
#' group label of the response variable is missing (i.e. NA)).
#' @param crit an integer indicating the impurity function used (1=Gini index / 2=Entropie/ 3=Misclassification rate).
#' @param penalty a boolean indicating if the decrease in node impurity must take account of the group size. Four penalty are
#' available: "No","Size","Root.size" or "Log".
#' @param case_min an integer indicating the minimun number of cases/non cases in a terminal nodes. The default is 1.
#' @param maxdepth the max depth for a split-tree.
#' @param IMPORTANCE a boolean indicating if the importance of each group need to be computed.
#'
#' @return a list with elements
#'    - tree : a data frame which summarizes the resulted CARTGV tree.
#'    - carts : a list containing all the CART objects used to buid the CARTGV tree.
#'    - splits : a list containing informations about the splits. Each element is an object retuned by the function "\code{split_cartgv}".
#'    - pop : a list containing the indices (rownames) of the observations which belong to the nodes.
#'    - tables_coupures : a list containing data frames that summarizes the splits.
#'    - importance_rand : a list providing the importance of each group at each node. Calculation of the importance is based on the Group Rand Importance.
#'    - importance_sur :  a list providing the importance of each group at each node. Calculation of the importance is based on the Group Surrogate Importance.
#'    - agreement : a list containing the measure of agreement between the selected split and the other possible split at each node and fr each group.
#'                This proximity measure is based on the Rand Index.
#'
#' @export
#'
cartgv<-function(data,
                 group,
                 crit=1,
                 case_min=1,
                 maxdepth=2,
                 penalty="No",
                 IMPORTANCE=TRUE)
  {

  ##Initialisation
  tree<-NULL
  parent<-c()# a vector indicatin g the ancestors of each node
  depth<-c()# a vector vector indicating th depth of eahc node
  var<-c() # a vector containing the indices of the splitting groups
  carts<-list() # a list containing all the CART objects used to buid the CARTGV tree
  n<-c() # a vector containing the number of observations in each node
  n_case<-c()# a vector containing the number of observations being "Y=1" in each node
  n_noncase<-c() # a vector containing the number of observations being "Y=0" in each node
  yval<-c() # a vector containing the node label = majority vote in each node
  prob<-c() #  a vector containing the prability  P[Y=1|N] in each node
  pop<-list() # a list containing the indices of the observations which belong to the nodes
  splits<-list() #a list containing informations about the splits
  tables_coupures<-list() #a list containing data frames that summarizes the splits
  improvment<-c() #Decrease in node impurity
  action<-c() # 1: if splitting; -1: at least one stopping criterion is fulfilled; -2: no impurity reduction
  igroups<-unique(group[!is.na(group)])
  if(IMPORTANCE==T){
    importance<-list()
    agreement<-list()
    importance_rand<-list()
    importance_sur<-list()
  }
  parent[1]<-NA
  depth[1]<-0
  groups_selec<-NULL
  n[1]<-dim(data)[1]
  pop[[1]]<-rownames(data)
  i<-1 #node index in the CARTGV tree. Remember that a node in CARTGV is a small CART tree named a splitting tree
  i_node_coupure<-NA #index of a node in the splitting tree
  yval[1] <- ifelse((length(which(data$Y == 1)) / n[1]) < 0.5, 0, 1)
  while(i<=length(n)){
    node<-data[intersect(unlist(pop[[i]]),rownames(data)),]
    prob[i]<-round((length(which(node$Y==1))/n[i]),4)
    n_case[i]<-length(which(node$Y==1))
    n_noncase[i]<-length(which(node$Y==0))
    if ((n_case[i] > case_min) & (n_noncase[i] > case_min)) {### one of the stopping criteria: node homogeneous?
        splits[[i]]<-split_cartgv(node=node,
                                group=group,
                                label=yval[i],
                                maxdepth=maxdepth,
                                penalty=penalty)

      improvment[i]<-max(unlist(splits[[i]][[crit]]))
      if (improvment[i] > 0){#Is there any group that reduces the node impurity
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
                                                                  unlist(splits[[i]]$carts[[ind_var]]$where)),
                                                            match.names=FALSE)$rand)

          importance_rand[[i]]<-unlist(importance[[i]])*unlist(agreement[[i]])

          imp_sur<-rep(NA,length(igroups))
          Ystar<-as.factor(splits[[i]]$carts[[ind_var]]$where)
          imp_sur[-ind_var]<-sapply(igroups[-ind_var],
                                    function(x)unlist(surrogate_split(Ystar,
                                                                      node,
                                                                      group=group,
                                                                      igroup=x,
                                                                      penalty=penalty,
                                                                      label=yval[i],
                                                                      maxdepth=maxdepth)[[crit]]))
          imp_sur[ind_var]<-max.importance
          importance_sur[[i]]<-imp_sur
        }
        var[i] <- igroups[ind_var]#Selection of the best split i.e. the group/split that maximizes the decrease in node impurity
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
        action[i]<--2 #No group/split reduces the impurity
        var[i]<-NA
        i<-i+1
      }
    }else{
      action[i]<--1 #There is at least one stopping criterion which is satisfied
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
  if(IMPORTANCE==FALSE){
    importance_rand<-NULL
    importance_sur<-NULL
    agreement<-NULL

  }
  return(list(tree=tree,
              carts=carts,
              splits=splits,
              pop=pop,
              tables_coupures=tables_coupures,
              importance_rand=importance_rand,
              importance_sur=cimportance_sur,
              agreement=agreement))
}

#' split_cartgv()
#'
#' Calculate the best split of a node for each group of input variables when building a CARTGV tree.
#'
#' @param node a data frame containing the observations in the node. The first column is the response vector named "Y" and with the lable
#' "0" and "1". The p-1 others variables are continuous or categroical variable must be coded as a set of dummy variables.
#' @param group a vector with the group number of each variable.
#' (WARNING : if there are "\code{p}" goups, the groups must be numbers from "\code{1}" to "\code{p}" in increasing order. The group label
#' of the response variable is missing (i.e. NA)).
#' @param label an integer indicating the label of the node (the majority class)
#' @param maxdepth an integer indicating the maximal depth for a split-tree. The default value is 2.
#' @param penalty a boolean indicating if the decrease in node impurity must take account of the group size. Four penalty are available: "No"
#' ,"Size","Root.size" or "Log".
#'
#' @return a list with elements
#'         - Gain_Gini: a vector containing the reduction of Gini in the node from splitting on each group,
#'         - Gain_Ent: a vector containing the reduction of Entropy in the node from splitting on each group,
#'         - Gain_Mis: a vector containing the reduction of the number of miclassified observations in the node from splitting on each group,
#'         - carts: a list containing for each group the CART object which summarizes the splitting tree,
#'         - pred: a matrix with "\code{nrows(node)}" lines and "\code{length(unique(group))}" columns containing for each group the prediction,
#'         resulting from the splitting tree.
#' @export
#'
#' @importFrom stats predict as.formula
split_cartgv<-function(node,
                       group,
                       label,
                       maxdepth=2,
                       penalty="No")
  {
  igroups<-unique(group[!is.na(group)])
  nb_group<-length(unique(group[which(group%in%igroups)]))
  Gain_Gini<-rep(0,nb_group)
  Gain_Ent<-rep(0,nb_group)
  Gain_Clas<-rep(0,nb_group)
  carts <- list()
  nb_nodes<-rep(0,nb_group)
  pred<-matrix(rep(NA,dim(node)[1]*nb_group),ncol=nb_group)
  for(j in 1:nb_group){
    ivar<-which(group==igroups[j])
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

    pred[,j]<-as.numeric.factor(stats::predict(cart, type = "class", newdata=node))
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
              pred=pred))
}


#' predict_cartgv
#'
#' Prediction of test data from a fitted CARTGV tree.
#'
#' The function called the function "\code{pred_cart}".
#'
#' @param new  a new data frame containing the same variables that "\code{data}".
#' @param tree the data frame "\code{tree}" returned by the function "\code{cartgv}".
#' @param carts the list "\code{carts}" returned by the function "\code{cartgv}".
#' @param coups the list "\code{table_coupures}" returned by the function "\code{cartgv}".
#'
#' @return a matrix with "\code{nrows(new)}" rows (the i-th row is provided prediction information about the -ith observation in "\code{new}")
#'         and the 4 followingcolumns,
#'         - hat.Y: the predicted label ("0" or "1"),
#'         - Y: the true label ("0" or "1"),
#'         - noeuds: name of the node in the CARTGV tree,
#'         - score:  conditionnal probability that the observation has "Y=1",
#'         - i_noeuds: name of the node in the splitting tree.
#'
#' @export
predict_cartgv<-function(new,tree, carts,coups){
  indx <- colnames(tree)
  indx <- indx[which(indx != 'leave')]
  tree[indx] <- lapply(tree[indx], as.numeric.factor)
  indx <- colnames(new)
  new[indx]  <- lapply(new[indx],  as.numeric.factor)

  predgv <- function( ind ) {
    i<-1
    while(tree$action[which(tree$node==i)]>=1){
      temp<-as.data.frame(t(pred_cart(ind,coups[[i]])[c(2,3,4,5)]))
      names(temp)<-c("i_node_coupure","pred","n_noncase","n_case")
      i<-tree$node[which( tree$parent==i
                          & tree$yval==temp$pred
                          & tree$n_case==temp$n_case
                          & tree$n_noncase==temp$n_noncase
                          & tree$i_node_coupure== temp$i_node_coupure)]
    }
    if(tree$action[which(tree$node==i)]<1){
      #on teste si le noeuds dans lequel tombre l'obs est une feuille
      return (c(tree$yval[which(tree$node==i)],            # pred (hat.Y)
               as.character(i),                            # noeuds
               tree$prob[which(tree$node==i)],             # score
               tree$i_node_coupure[which(tree$node==i)]))  # i_noeuds

      i<-1
    } else {
      return(c(NA,NA,NA,NA))
  }
  }

  P<-dim(new)[1]

  res <- matrix()
  length(res) <- 5 * P
  dim(res) <- c(P, 5)
  res[,1] <- new[,1]
  res[,2:5] <- t(apply(new,1,predgv))

  colnames(res)<-c("Y","hat.Y","noeuds","score","i_noeuds")
  if (2 %in% res[,1]) {res[,1] <- ifelse(res[,1]==2,1,0)}
  if (2 %in% res[,2]) {res[,2] <- ifelse(res[,2]==2,1,0)}
  return(as.data.frame(res))
}




#' impurity.cartgv
#'
#' Compute the impurity information from a sequence of subtrees by using an independent set.
#' This function is used to select the best subtree among a sequence of subtrees.
#'
#' @param validation a new data frame containing the same variables that "\code{data}".
#' @param tree_seq the object returned by the function "\code{extract_subtrees}".
#' It is a sequence of optimal subtrees obtained by applying
#' the cost-complexity pruning method
#' @param tree a fitted CARTGV tree. It is an object returned by the function "\code{cartgv}".
#'
#' @return a list with elements
#'    - impurete: a data frame containing the value of several impurity fucntions
#'        (in this order Gini, Entropy, misclassification rate)
#'                    for each subtree of the sequence. The i-th row corresponds to the i-th
#'                    subtree of the sequence.
#'    - pred:  a list containing the prediction of the label for the data set "\code{validation}"
#'       based on each subtree. Precisely, the i-th element is the object returned by the function
#'                 "\code{predict_cartgv}" for the i-th subtree nd by using
#'                 the data set "\code{validation}".
#'    - summary_noeuds: a list containg for each subtree informations about the nodes
#'        (nom_noeuds: node name, N: number of observations in
#'                          the node, \code{N[Y=1]}: number of observation with "Y=1" in
#'                          the node, \code{N[Y=0]}: number of observation with "Y=0"
#'                          in the node, \code{P[Y=1]}: estimated probability that an
#'                          observation in the node is assigned to the label "Y=1",
#'                          \code{P[Y=1]}: estimated probability that an observation in
#'                          the node is assigned to the label "Y=0" and \code{P[hat.Y!=Y]}:
#'                          misclassification rate in the node).
#' @export
#'
impurity.cartgv <- function(validation, tree_seq,tree) {
  N_tree <- length(tree_seq)
  N <- dim(validation)[1]
  pred <- list()
  sum_noeuds <- list()
  impurete <- matrix(rep(NA, N_tree * 3), nrow = N_tree, ncol = 3)
  for (k in 1:N_tree) {
    predictions <-as.data.frame(predict_cartgv(validation,
                                               as.data.frame(tree_seq[[k]]),
                                               tree$carts,
                                               tree$tables_coupures))

    predictions$noeuds  <- as.numeric.factor(predictions$noeuds)
    predictions$Y       <- as.numeric.factor(predictions$Y)

    predictions<-predictions[order(predictions$noeuds),]
    n <- as.numeric(table(predictions$noeuds))
    n1 <- tapply(predictions$Y, predictions$noeuds,sum)
    p1 <- n1 / n
    p0 <- 1 - p1
    predictions$error.pred <-as.numeric.factor(apply(predictions[, c("hat.Y", "Y")],
                                                           1,
                                                           function(x)ifelse(x[1] != x[2], 1, 0)))
    misclass <-tapply(predictions$error.pred,
                      predictions$noeuds,
                      sum)

    summaryNoeuds <-as.data.frame(cbind(as.numeric.factor(names(table(predictions$noeuds))),
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

  list(impurete = impurete,
       pred = pred,
       summary_noeuds = sum_noeuds)
}


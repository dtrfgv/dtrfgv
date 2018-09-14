############################################################################################
###########################        Tools for rpart objects       ########################### 
############################################################################################ 

as.numeric.factor <- function(x) {as.numeric(as.character(x))}

#' calcul_cart
#'
#' Create a data frame that sums up a CART tree obtained with the function in the rpart package
#'
#' WARNING: the option 'control=rpart.control(maxsurrogate=0, maxcompte=0)' must be used in the function rpart when building the rpart object. 
#'           Otherwise the result returned by the function may be incorrect
#'
#' @param cart  an rpart object
#' @param data  the data used to create the rpart object. It must be a data frame
#'
#' @return tableau  a data frame that sums up the CARt tree:
#'                  - node_name_cart: node number
#'                  - parent: ancestor of the node 
#'                  - depth: depth of the node
#'                  - var: name of the splitting variables
#'                  - threshold: splitting value 
#'                  - n: number of obsevations in the node 
#'                  - pred: label of the node
#'                  - n_case: number of observation with "\code{Y=1}" in the node 
#'                  - n_noncase: number of observation with "\code{Y=0}" in the node 
#'                  - \code{prob = P[Y=1|N]}: empirical probability that in the node an observation belongs to the label "\code{Y=1}" 
#'                  - sens: direction for the decision rule. Observation with "\code{x<threshold}" are sent to the left node and observations with "\code{x>=threshold}" are sent to the right node.
#' 
#' @import rpart
#'
calcul_cart <- function(cart,data) {
  table<-cart$frame
  split<-cart$splits
  if(cart$control$maxcompete!=0 & cart$control$maxsurrogate!=0){
    print("error: you need to put control=rpart.control(maxsurrogate=0, maxcompete=0) in rpart()")
    tableau<-NULL
  }else{
    if(dim(table)[1]>1){
      N<-dim(table)[1]
      node_cart<-rownames(table)
      max<-max(as.numeric.factor(node_cart))
      depth<-c(0,rep(NA,N-1))
      parent<-rep(NA,N)
      threshold<-rep(NA,N)
      sens<-rep(NA,N)
      d<-1
      while((2^d)<=max){
        i_nodes<-which(  (as.numeric.factor(node_cart) <  2^(d+1)) 
                       & (as.numeric.factor(node_cart) >= 2^(d) ))
        parents<-rep(seq(2^(d-1),2^(d)-1,1),rep(2,2^(d)-(2^(d-1))))
        nodes<-seq(2^(d),2^(d+1)-1,1)
        depth[i_nodes]<-d
        parent[i_nodes]<-sapply(
          i_nodes,function(x){parents[which(nodes==as.numeric.factor(node_cart[x]))]})
        d<-d+1
      }
      tableau<-as.data.frame(cbind(node_cart,
                                   parent,
                                   depth,
                                   as.character(table$var),
                                   threshold,table$n,
                                   as.numeric.factor(table$yval)-1,
                                   table$yval2[,c(2,3,5)]))
      
      names(tableau)<-c("node_name_cart",
                        "parent",
                        "depth",
                        "var",
                        "threshold",
                        "n",
                        "pred",
                        "n_noncase",
                        "n_case",
                        "prob")
      
      tableau$threshold<-as.numeric.factor(tableau$threshold)
      tableau$threshold[which(tableau$var!="<leaf>")]<-as.numeric.factor(split[,"index"])
      sens[which(tableau$var!="<leaf>")]<-ifelse(as.numeric.factor(split[,"ncat"])==-1,"inf","sup")
      tableau$sens<-sens
      tableau
      
    }else{
      print("no split: the final tree is trivial")
      tableau<-as.data.frame(cbind(1,
                                   NA,
                                   0,
                                   "<leaf>",
                                   NA,
                                   table$n,
                                   as.numeric.factor(table$yval)-1,
                                   table$yval2[2],
                                   table$yval2[3],
                                   table$yval2[5],
                                   NA))
      
      names(tableau)<-c("node_name_cart",
                        "parent",
                        "depth",
                        "var",
                        "threshold",
                        "n",
                        "pred",
                        "n_noncase",
                        "n_case",
                        "prob",
                        "sens")
    }
  }
  return(tableau)
}


#' pred_cart
#' 
#' Function used to predict the label ("0" or "1") of a new observation from a fitted CARTGV object. 
#' 
#' @param new_obs a new observation. It must be a vector containing the values of the variables of the learning data set used to built the rpart object
#' @param tableau object returns by the function 'calcul_cart' 
#'
#' @return a vector with elements (in this order): number of the predicted node in the splitting tree, index of the predicted node in the object returned 
#'                                                 by the function 'pred_cart', the label of the predicted node, the number of observations with "\code{Y=1}" 
#'                                                 in the predicted node, the number of observations with "\code{Y=0}" in the predicted node and the empirical 
#'                                                 probability that in the predicted node an observation belongs to the label "\code{Y=1}" 
#'
#'
#'
pred_cart<-function(new_obs,tableau){
  i_node<-1
  while(as.character(tableau$var[i_node])!="<leaf>"){
    i_sons<-which(as.numeric.factor(tableau$parent)==as.numeric.factor(tableau$node_name_cart[i_node]))
    if(new_obs[as.character(tableau$var[i_node])]>=as.numeric.factor(tableau$threshold[i_node])){
      if(tableau$sens[i_node]=="sup"){
        i_node<-i_sons[1]
      }else{
        i_node<-i_sons[2]
      }
    }else{    
      if(tableau$sens[i_node]=="sup"){
        i_node<-i_sons[2]
      }else{
        i_node<-i_sons[1]
      }
    }
  }
  pred_node_name<-as.character(tableau$node_name_cart[i_node])
  pred_i_node<-as.numeric.factor(i_node)
  return(c(as.numeric.factor(pred_node_name),
           as.numeric.factor(pred_i_node),
           as.numeric.factor(tableau$pred[i_node]),
           as.numeric.factor(tableau$n_noncase[i_node]),
           as.numeric.factor(tableau$n_case[i_node]),
           as.numeric.factor(tableau$prob[i_node])))
}




#' impurete_rpart
#' 
#' Compute the impurity of a sequence of trees which are based on a rpart object.  
#' 
#' a partir d'un echantillon independant et d'une sequence d'arbres emboites
#' La fonction prend en entree un échantillon test et une sequence d'arbres emboites
#' Pour chaque arbre, on predit la classe de chaque observation de l'ensemble test;
#' puis a partir de ces resultats, on calcul l'impurete de chaque arbre (gini, l'entropie et 
#' le taux de mal-classes sont calcules)
#'
#' @param validation an new data set. It must be a data frame containing the same variables that those contained in the learning data set used to built the rpart object.
#' @param tree_seq   a sequence of subtrees. It must be a list where each element is an object returned by the function 'rpart::snip.rpart' 
#'
#' @return a list with elements
#'  - impurete:  a matrix containing the impurity values (respectively Gini, Entropy and Misclassification rate) of each subtrees evaluated on the data set '\code{validation}'.
#'  - pred:      a list containing the prediction vector of each subtree
#'  - summary_noeuds:  a list providing informations about each subtree:
#'      - nom_noeuds: number of the node
#'      - N: number of observations in the node
#'      - \code{N[Y=1]}: number of observation with "\code{Y=1}" in the node 
#'      - \code{P[Y=1]}: empirical probability that in the node an observation belongs to the label "\code{Y=1}" 
#'      - \code{P[Y=0]}: empirical probability that in the node an observation belongs to the label "\code{Y=0}"
#'      - \code{P[hat.Y!=Y]}: misclassification rate in the node
#'
impurete_rpart <- function(validation, tree_seq) {
  N_tree <- length(tree_seq)
  N <- dim(validation)[1]
  pred <- list()
  sum_noeuds <- list()
  
  validation$Y <- as.numeric.factor(validation$Y)
  
  impurete <- matrix(rep(NA, N_tree * 3), nrow = N_tree, ncol = 3)
  for (k in 1:N_tree) {
    p <- stats::predict(tree_seq[[k]], type = "matrix", validation)
    noeuds <-apply(p, 1, 
                   function(x)rownames(tree_seq[[k]]$frame[which( (tree_seq[[k]]$frame$yval2[, 6] == x[6]) 
                                                                 &(tree_seq[[k]]$frame$yval2[, 2] == x[2]) 
                                                                 &(tree_seq[[k]]$frame$yval2[, 3] == x[3]) 
                                                                 &(tree_seq[[k]]$frame$yval2[, 4] == x[4]) 
                                                                 &(tree_seq[[k]]$frame$yval2[, 5] == x[5]) 
                                                                 &(tree_seq[[k]]$frame$yval==x[1]))[1], ]))
    
    predictions <-as.data.frame(cbind(p[, 1] - 1, validation$Y, noeuds))
    names(predictions) <- c("hat.Y", "Y", "noeuds")
    n <- as.numeric(table(as.numeric.factor(predictions$noeuds)))
    n1 <- tapply(as.numeric.factor(predictions$Y), 
                 as.numeric.factor(predictions$noeuds),sum)
    p1 <- n1 / n
    p0 <- 1 - p1
    predictions$error.pred <-as.numeric.factor(apply(predictions[, c("hat.Y", "Y")], 
                                                           1, 
                                                           function(x)ifelse(x[1] != x[2], "1", "0")))
    misclass <-tapply(as.numeric.factor(predictions$error.pred), 
                      as.numeric.factor(predictions$noeuds), sum)
    summaryNoeuds <-as.data.frame(cbind(as.numeric.factor(names(table(as.numeric.factor(predictions$noeuds)))), 
                                        n, n1, p1, p0, misclass))
    names(summaryNoeuds) <-c("nom_noeuds","N","N[Y=1]","P[Y=1]","P[Y=0]","P[hat.Y!=Y]")
    impurete[k, 1] <- sum(apply(summaryNoeuds, 1, function(x)x[2] * gini(x[c(4, 5)]))) / N
    impurete[k, 2] <- sum(apply(summaryNoeuds, 1, function(x)x[2] * entropy(x[c(4, 5)]))) / N
    impurete[k, 3] <- sum(summaryNoeuds[, 6]) / (N)
    impurete<-as.data.frame(impurete)
    names(impurete)<-c("Gini","Information","Misclass")
    pred[[k]]<-predictions
    sum_noeuds[[k]]<-summaryNoeuds
    predictions<-NULL
    summaryNoeuds<-NULL
  }
  list(impurete = impurete,pred = pred,summary_noeuds = sum_noeuds)
}

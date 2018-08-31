######## Manipulation des arbres CART (package rpart)

library(rpart)

##################################################################################################################################################################
# Reconstruction  des arbres CART
##################################################################################################################################################################

#La fonction construit un data.frame qui permet de donner toutes les infos necessaire à la reconstruction de l'arbre cart construit par RPART
### Paramètres d'entrée :
# cart : objet rpart
# data  : data.frame ayant servi à la construction de l'objet cart
### Attention cette fonction n'est valide que si rpart est parametré 
### avec control=rpart.control(maxsurrogate=0, maxcompte=0)
#' Title
#'
#' @param cart 
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
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
      max<-max(as.numeric(as.character(node_cart)))
      depth<-c(0,rep(NA,N-1))
      parent<-rep(NA,N)
      threshold<-rep(NA,N)
      sens<-rep(NA,N)
      d<-1
      while((2^d)<=max){
        i_nodes<-which((as.numeric(as.character(node_cart))<(2^(d+1))) 
                       & (as.numeric(as.character(node_cart))>=(2^(d))))
        parents<-rep(seq(2^(d-1),2^(d)-1,1),rep(2,2^(d)-(2^(d-1))))
        nodes<-seq(2^(d),2^(d+1)-1,1)
        depth[i_nodes]<-d
        parent[i_nodes]<-sapply(
          i_nodes,function(x){parents[which(nodes==as.numeric(as.character(node_cart[x])))]})
        d<-d+1
      }
      tableau<-as.data.frame(cbind(node_cart,parent,depth,as.character(table$var),threshold,table$n,
                                   as.numeric(as.character(table$yval))-1,table$yval2[,c(2,3,5)]))
      names(tableau)<-c("node_name_cart","parent","depth","var","threshold","n","pred","n_noncase","n_case","prob")
      tableau$threshold<-as.numeric(as.character(tableau$threshold))
      tableau$threshold[which(tableau$var!="<leaf>")]<-as.numeric(as.character(split[,"index"]))
      sens[which(tableau$var!="<leaf>")]<-ifelse(as.numeric(as.character(split[,"ncat"]))==-1,"inf","sup")
      tableau$sens<-sens
      # i_leafs<-which(tableau$var=="<leaf>")
      # for(l in 1:length(i_leafs)){
      #   i_node<-i_leafs[l]
      #   ind<-data[which(cart$where==i_node),]#cart$where donne l'indice de la feuille et non le nom du noeud
      #   while(!is.na(tableau$parent[i_node])){
      #     i_parent<-which(as.numeric(as.character(tableau$node_name_cart))==
      #                       as.numeric(as.character(tableau$parent[i_node])))
      #     if(mean(ind[,as.character(tableau$var[i_parent])])>=
      #        as.numeric(as.character(tableau$threshold[i_parent]))){
      #       sens[i_node]<-"sup"
      #     }else{
      #       sens[i_node]<-"inf"
      #     }
      #     i_node<-i_parent
      #   }
      # }
      # tableau$sens <- sens
      tableau
      
    }else{
      print("no split: the final tree is trivial")
      tableau<-as.data.frame(cbind(1,NA,0,"<leaf>",NA,table$n,
                                   as.numeric(as.character(table$yval))-1,table$yval2[2],table$yval2[3],table$yval2[5],NA))
      
      names(tableau)<-c("node_name_cart","parent","depth","var","threshold","n","pred","n_noncase","n_case","prob","sens")
    }
  }
  return(tableau)
}
# Valeurs retourn'ees par la fonction : 
#==> la fonction renvoie une liste donnant : 
###   tableau : un data.frame qui donne toutes les infos sur l'arbre cart
###             ("node_name_cart"= numero du noeud,"parent"=ancêtre,"depth"=profondeur du noeud,
###              "var"=variable utilisée pour couper le noeud,"threshold"=seuil de la coupure,
###               "n"=effectif du noeud,"pred"=label du neoud,"n_noncase"= nombre d'observation "Y=0" dans le noeuds
###               ,"n_case"=nombre d'individus "Y=1" dans le noeud,"prob"= P[Y=1|N], 
###               "sens"=sens pour la règle de coupure divisant son noeud parent: 
###                     le sous ensembles des observations verifiant x<threshold est envoyé à gauche si "inf" à droite sinon. )



##################################################################################################################################################################
# Prediction de la feuille
##################################################################################################################################################################

### Foncion permettant de donner le nom cart et l'indice du noeud (sortie de la fonction carts$where)
### qui contient une nouvelle observation
### Paramètres d'entrée:
# new_obs : une nouvelle observation (avec les mêmes variables que celle de l'echantillonn d'apprentissage)
# tableau : sortie de la fonction calcul_cart()
pred_cart<-function(new_obs,tableau){
  i_node<-1
  while(as.character(tableau$var[i_node])!="<leaf>"){
    #print(paste("i_node:",i_node))
    i_sons<-which(as.numeric(as.character(tableau$parent))==as.numeric(as.character(tableau$node_name_cart[i_node])))
    #print(paste("i_sons:",i_sons))
    if(new_obs[as.character(tableau$var[i_node])]>=as.numeric(as.character(tableau$threshold[i_node]))){
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
  pred_i_node<-as.numeric(as.character(i_node))
  return(c(as.numeric(as.character(pred_node_name)),as.numeric(as.character(pred_i_node)),as.numeric(as.character(tableau$pred[i_node])),
           as.numeric(as.character(tableau$n_noncase[i_node])),as.numeric(as.character(tableau$n_case[i_node])),as.numeric(as.character(tableau$prob[i_node]))))
}



##################################################################################################################################################################
# Elagage :  calcul de l'impuret'e a' partir d'un 'echantillon independant et d'une s'equence d'arbres emboites
##################################################################################################################################################################

#La fonction prend en entrée un échantillon test et une séquence d'arbres emboités
#Pour chaque arbre, on prédit la classe de chaque observation de l'ensemble test;
#puis à partir de ces résultats, on calcul l'impureté de chaque arbre (gini, l'entropie et le taux de mal-classés sont calculés)

impurete_rpart <- function(validation, tree_seq) {
  N_tree <- length(tree_seq)
  N <- dim(validation)[1]
  pred <- list()
  sum_noeuds <- list()
  impurete <- matrix(rep(NA, N_tree * 3), nrow = N_tree, ncol = 3)
  for (k in 1:N_tree) {
    p <- predict(tree_seq[[k]], type = "matrix", validation)
    noeuds <-apply(p, 1, function(x)rownames(tree_seq[[k]]$frame[which((tree_seq[[k]]$frame$yval2[, 6] == x[6]) &(tree_seq[[k]]$frame$yval2[, 2] == x[2]) &(tree_seq[[k]]$frame$yval2[, 3] == x[3]) &
                                                                         (tree_seq[[k]]$frame$yval2[, 4] == x[4]) & (tree_seq[[k]]$frame$yval2[, 5] == x[5]) & (tree_seq[[k]]$frame$yval==x[1]))[1], ]))
    predictions <-as.data.frame(cbind(p[, 1] - 1, as.numeric(as.character(validation$Y)), noeuds))
    names(predictions) <- c("hat.Y", "Y", "noeuds")
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
    sum_noeuds[[k]]<-summaryNoeuds
    predictions<-NULL
    summaryNoeuds<-NULL
  }
  list(impurete = impurete,pred = pred,summary_noeuds = sum_noeuds)
}
# Valeurs retourn'ees par la fonction : 
#==> la fonction renvoie une liste donnant : 
###   impurete : une matrice contenant les differentes valeurs d'impuret'e pour chaque sous-arbres
###   pred = liste des predictions pour chaque sous arbres
###   summary_noeuds : liste contenant pour chaque sous-arbre des infos sur ses noeuds: 
###                     nom_noeuds = noms du noeuds
###                     N = nb obs dans le neoud
###                     N[Y=1] = nb d'obs dans le noeuds appartenant a' la class "Y=1" 
###                     P[Y=1] = proba pour une obs du noeuds d'appartenanir  a' la class "Y=1"
###                     P[Y=0] = proba pour une obs du noeuds d'appartenanir  a' la class "Y=0"
###                     P[hat.Y!=Y] = tx de mal class'es dans le noeud

# =========================================================================================
# predict.cartgv()
# =========================================================================================


### la fonction fait appel à la fonction pred_cart()
### IMPORTANT : i_noeuds (indice du noeud dans l'arbre de coupure, càd indice donné par rpart) et noeuds (indice du noeud dans l'arbre cartgv) utilisées comme clés primaires pour identifier de manière unique un noeuds
#' Title
#'
#' @param new 
#' @param tree 
#' @param carts 
#' @param coups 
#'
#' @return
#' @export
#'
#' @examples
predict.cartgv<-function(new,tree,carts,coups){
  P<-dim(new)[1]
  pred<-rep(NA,P)
  noeuds <- rep(NA, P)
  i_noeuds <- rep(NA, P)
  score <- rep(NA, P)
  for(p in 1:P){
    ind<-new[p,]
    i<-1
    while(as.numeric(as.character(tree$action[which(as.numeric(as.character(tree$node))==i)]))>=1){
      temp<-as.data.frame(t(pred_cart(ind,coups[[i]])[c(2,3,4,5)]))
      names(temp)<-c("i_node_coupure","pred","n_noncase","n_case")
      i<-as.numeric(as.character(tree$node[which(as.numeric(as.character(tree$parent))==i
                                                 & as.numeric(as.character(tree$yval))==(as.numeric(as.character(temp$pred)))
                                                 & as.numeric(as.character(tree$n_case))==as.numeric(as.character(temp$n_case))
                                                 & as.numeric(as.character(tree$n_noncase))==as.numeric(as.character(temp$n_noncase))
                                                 & as.numeric(as.character(tree$i_node_coupure))== as.numeric(as.character(temp$i_node_coupure)))]))
    }
    if(as.numeric(as.character(tree$action[which(as.numeric(as.character(tree$node))==i)]))<1){
      #on teste si le noeuds dans lequel tombre l'obs est une feuille 
      pred[p]<-as.numeric(as.character(tree$yval[which(as.numeric(as.character(tree$node))==i)]))
      noeuds[p] <- i
      i_noeuds[p]<-as.numeric(as.character(tree$i_node_coupure[which(as.numeric(as.character(tree$node))==i)]))
      score[p]<- as.numeric(as.character(tree$prob[which(as.numeric(as.character(tree$node))==i)]))
      i<-1
    }
  }
  res<-as.data.frame(cbind(Y=as.numeric(as.character(new$Y)),hat.Y=as.numeric(as.character(pred)),noeuds=as.character(noeuds),score=as.numeric(as.character(score)),i_noeuds=as.numeric(as.character(i_noeuds))))
  names(res)<-c("Y","hat.Y","noeuds","score","i_noeuds")
  if("2"%in%res$Y){res$Y<-ifelse(res$Y=="2","1","0")}
  if("2"%in%res$hat.Y){res$hat.Y<-ifelse(res$hat.Y=="2","1","0")}
  return(res)
}


# Valeurs retourn'ees par la fonction : 
#==> la fonction renvoie une matrice donnant pour chaque individu de new (dans l'ordre) : 
###   hat.Y :la classe predite
###   Y : la classe d'appartenance
###   noeuds : le numero du noeud contenant l'observation
###   score : le score le l'observation


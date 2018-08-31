##################################################################################################################################################################
# Elagage :  calcul de l'impuret'e a' partir d'un 'echantillon independant et d'une s'equence d'arbres emboites

# WARNING : fonction corrigée le 19 janvier2018 car il y avait un problème d'ordre pour les noeuds dans les apply
##################################################################################################################################################################

#La fonction prend en entrée un échantillon test et une séquence d'arbres emboités
#Pour chaque arbre, on prédit la classe de chaque observation de l'ensemble test;
#puis à partir de ces résultats, on calcul l'impureté de chaque arbre (gini, l'entropie et le taux de mal-classés sont calculés)
### Paramètres d'entrée :
# validation : data.frame contenant les mêmes variables que le jeux de données utilisés pour construire l'objet tree
# tree  : objet retourné par la fonction Tree_CART_bin() ou Tree_CART()
# method : la nature des coupures pour l'arbre CC binaire ("bin" si tree construit avec Tree_CART_bin()) 
#          ou ("multi" si tree construit avec Tree_CART()) 

#' Title
#'
#' @param validation 
#' @param tree_seq 
#' @param tree 
#'
#' @return
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
    predictions <-as.data.frame(predict.cartgv(validation,as.data.frame(tree_seq[[k]]),tree$carts,tree$tables_coupures))
    predictions<-predictions[order(as.numeric(as.character(predictions$noeuds))),]
    n <- as.numeric(table(as.numeric(as.character(predictions$noeuds))))
    n1 <- tapply(as.numeric(as.character(predictions$Y)), as.numeric(as.character(predictions$noeuds)),sum)
    p1 <- n1 / n
    p0 <- 1 - p1
    predictions$error.pred <-as.numeric(as.character(apply(predictions[, c("hat.Y", "Y")], 1, function(x)ifelse(x[1] != x[2], "1", "0"))))
    misclass <-tapply(as.numeric(as.character(predictions$error.pred)), as.numeric(as.character(predictions$noeuds)), sum)
    summaryNoeuds <-as.data.frame(cbind(as.numeric(as.character(names(table(as.numeric(as.character(predictions$noeuds)))))), n, n1, p1, p0, misclass))
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


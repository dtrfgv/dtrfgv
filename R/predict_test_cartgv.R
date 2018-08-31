# =========================================================================================
# predict.test.cartgv() ### Test Pierre
# =========================================================================================

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

# Valeurs retourn'ees par la fonction :
#==> la fonction renvoie une matrice donnant pour chaque 
#==> individu de new (dans l'ordre) :
###   - hat.Y :la classe predite
###   - Y : la classe d'appartenance
###   - noeuds : le numero du noeud contenant l'observation
###   - score : le score le l'observation



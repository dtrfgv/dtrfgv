###############################################################################
# Prediction de la feuille
##############################################################################

# Fonction permettant de donner le nom cart et l'indice du noeud (sortie de la 
# fonction carts$where)
# qui contient une nouvelle observation
# Paramètres d'entrée:
# new_obs : une nouvelle observation (avec les mêmes variables que celle de 
#           l'echantillon d'apprentissage)
# tableau : sortie de la fonction calcul_cart()
pred_cart <- function(new_obs, tableau) {
  i_node <- 1
  while (as.character(tableau$var[i_node]) != "<leaf>") {
    #print(paste("i_node:",i_node))
    i_sons <-
      which(as.numeric(as.character(tableau$parent)) == as.numeric(as.character(tableau$node_name_cart[i_node])))
    if (new_obs[as.character(tableau$var[i_node])] >= as.numeric(as.character(tableau$threshold[i_node]))) {
      if (tableau$sens[i_node] == "sup") {
        i_node <- i_sons[1]
      } else{
        i_node <- i_sons[2]
      }
    } else{
      if (tableau$sens[i_node] == "sup") {
        i_node <- i_sons[2]
      } else{
        i_node <- i_sons[1]
      }
    }
  }
  pred_node_name <- as.character(tableau$node_name_cart[i_node])
  pred_i_node <- as.numeric(as.character(i_node))
  return(c(
    pred_node_name,
    pred_i_node,
    as.numeric(as.character(tableau$pred[i_node])),
    as.numeric(as.character(tableau$n_noncase[i_node])),
    as.numeric(as.character(tableau$n_case[i_node])),
    as.numeric(as.character(tableau$prob[i_node]))
  ))
}

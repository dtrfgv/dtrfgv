#' perf
#' 
#' Optimal cut
#' 
#' Trouver le seuil qui permet le meilleur compromis sensibilite-specificite
#'
#' @param perf 
#'
#' @return NULL
#'
opt.cut = function(perf) {
  cut.ind = mapply(
    FUN = function(x, y, p) {
      d = (x - 0) ^ 2 + (y - 1) ^ 2
      ind = which(d == min(d))
      if (length(ind) == 1) {
        c( sensitivity = y[[ind]],
           specificity = 1 - x[[ind]],
           cutoff = p[[ind]])
      }else{
        c(sensitivity = y[[ind[length(ind)]]],
          specificity = 1 - x[[ind[length(ind)]]],
          cutoff = p[[ind[length(ind)]]])
      }
    }
    ,
    perf@x.values,
    perf@y.values,
    perf@alpha.values
  )
}


#' performances
#'
#' @param score 
#' @param test 
#'
#' @return
#' @importFrom ROCR performance
performances <- function(score, test) {
  pred <- ROCR::prediction(score, as.factor(test[, 1]))
  ## ROC curve
  roc <- performance(pred, measure = "tpr", x.measure = "fpr")
  ## Best cut-off
  c <- opt.cut(roc)
  ##
  predY <- ifelse(score >= 0, "1", "0")
  if (length(unique(predY)) == 1) {
    if (unique(predY) == "0") {
      m <- table(predY, as.factor(test[, 1]))
      tab <- as.table(matrix(c(m[1, 1], m[1, 2], 0, 0), byrow = T, nrow =2))
      rownames(tab) <- c("0", "1")
      colnames(tab) <- c("0", "1")
      xtab <- confusionMatrix(tab, positive = "1")
    } else
    {
      m <- table(predY, as.factor(test[, 1]))
      tab <- as.table(matrix(c(0, 0, m[1, 1], m[1, 2]), byrow = T, nrow =2))
      rownames(tab) <- c("0", "1")
      colnames(tab) <- c("0", "1")
      xtab <- confusionMatrix(tab, positive = "1")
    }
  } else{
    xtab <- confusionMatrix(table(predY, as.factor(test[, 1])), positive = "1")
  }
  table <- xtab$table
  conc  <- (sum(diag(table(predY, as.factor(test[, 1])))) / sum(table(predY, as.factor(test[, 1]))))
  sens  <- as.numeric(xtab$byClass[1])
  spe   <- as.numeric(xtab$byClass[2])
  auc   <- performance(pred, "auc")@y.values[[1]]
  
  resultats <- list(table, conc, sens, spe, auc, c[3])
  return(resultats)
}



#' xtab function
#' 
#' Matrice de confusion + performances predictives 
#' 
#' (gere quand le modele predit tout en 1 ou 0)
#'
#' @param predY 
#' @param Y 
#'
#' @return
#' 
#' @importFrom caret confusionMatrix
#'
xtab_function<-function(predY,Y){
  if (length(unique(predY)) == 1) {
    if (unique(predY) == "0") {
      m <- table(predY, as.factor(Y))
      tab <- as.table(matrix(c(m[1, 1], m[1, 2], 0, 0), byrow = T, nrow = 2))
      rownames(tab) <- c("0", "1")
      colnames(tab) <- c("0", "1")
      xtab <- confusionMatrix(tab, positive = "1")
    } else
    {
      m <- table(predY, as.factor(Y))
      tab <- as.table(matrix(c(0, 0, m[1, 1], m[1, 2]), byrow = T, nrow = 2))
      rownames(tab) <- c("0", "1")
      colnames(tab) <- c("0", "1")
      xtab <- confusionMatrix(tab, positive = "1")
    }
  } else{
    xtab <- confusionMatrix(table(as.factor(predY), as.factor(Y)), positive = "1")
  }
  xtab
}

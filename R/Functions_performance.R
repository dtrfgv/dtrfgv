############################################################################################
########     Functions to calculate the performance of a binary classifie  #################
############################################################################################


#' opt.cut
#' 
#' Find optimal cutoff (for the trade-off Sensitivity/specificity)
#'
#' @param perf an S4 object of class performance with elements to create the ROC curve  (y.values: True positive rate, x.values: False positive rate and alpha.value: the threshold)
#'
#' @return 
#' cut.ind : the index of the optimal cut
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
#' Calculates a cross-tabulation of observed and predicted classes with associated statistics for the optimal cutoff (i.e. the cutoff that makes the best trade-off between sensitivity and specificity) 
#'
#' As the function xtab_function, this function deals with the case where a model assigns all observation to the same class (0 or 1)
#'
#' @param score  probabilistic output 
#' @param Y      a vector containing the true class labels. Must have the same dimensions as ’score’. It must be a factor with levels "0" and "1"
#'
#' @return resultats a list with elements:
#'                   - table:   the confusion matrix (i.e. a cross-tabulation of observed and predicted classes) for the optimal cutoff
#'                   - conc:    the accuracy rate for the optimal cutoff
#'                   - sens:    the sensitivy rate for the optimal cutoff
#'                   - spe:     the specificity rate for the optimal cutoff
#'                   - auc:     the auc for the optimal cutoff
#'                   - cutoff:  the optimal cutoff
#'                    
#' @importFrom ROCR performance
performances <- function(score, Y) {
  pred <- ROCR::prediction(score, Y)
  ## ROC curve
  roc <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
  ## Best cut-off
  cutoff <- opt.cut(roc)[3]
  ##
  predY <- ifelse(score >= cutoff, "1", "0")
  if (length(unique(predY)) == 1) {
    if (unique(predY) == "0") {
      m <- table(predY, Y)
      tab <- as.table(matrix(c(m[1, 1], m[1, 2], 0, 0), byrow = T, nrow =2))
      rownames(tab) <- c("0", "1")
      colnames(tab) <- c("0", "1")
      xtab <- caret::confusionMatrix(tab, positive = "1")
    } else
    {
      m <- table(predY, Y)
      tab <- as.table(matrix(c(0, 0, m[1, 1], m[1, 2]), byrow = T, nrow =2))
      rownames(tab) <- c("0", "1")
      colnames(tab) <- c("0", "1")
      xtab <- caret::confusionMatrix(tab, positive = "1")
    }
  } else{
    xtab <- caret::confusionMatrix(table(predY, Y), positive = "1")
  }
  table <- xtab$table
  conc  <- xtab$overall[1]
  sens  <- as.numeric(xtab$byClass[1])
  spe   <- as.numeric(xtab$byClass[2])
  auc   <- performance(pred, "auc")@y.values[[1]]
  
  resultats <- list(table, conc, sens, spe, auc, cutoff)
  return(resultats)
}



#' xtab function
#' 
#' Calculates a cross-tabulation of observed and predicted classes with associated statistics
#'
#' Compared to the function caret::confusionMatrix, this function deals with the case where a model assigns all observation to the same class (0 or 1)
#' 
#'
#' @param predY A vector containing the predictions
#' @param Y     A vector containing the true class labels. Must have the same dimensions as ’predY’
#'
#' @return xtab A list with elements
#'    xtab A list with elements: , , , 
#'    - table : the results of table on data and reference
#'    - positive : the positive result level
#'    - overall : a numeric vector with overall accuracy and Kappa statistic values
#'    - byClass : the sensitivity, specificity, positive predictive value, negative predictive value, precision, recall, F1, 
#'                prevalence, detection rate, detection prevalence and balanced accuracy for each class. For two class systems, 
#'                this is calculated once using the positive argument 
#' @export

#' @importFrom caret confusionMatrix
#'
xtab_function<-function(predY,Y){
  if (length(unique(predY)) == 1) {
    if (unique(predY) == "0") {
      m <- table(predY, as.factor(Y))
      tab <- as.table(matrix(c(m[1, 1], m[1, 2], 0, 0), byrow = T, nrow = 2))
      rownames(tab) <- c("0", "1")
      colnames(tab) <- c("0", "1")
      xtab <- caret::confusionMatrix(tab, positive = "1")
    } else
    {
      m <- table(predY, as.factor(Y))
      tab <- as.table(matrix(c(0, 0, m[1, 1], m[1, 2]), byrow = T, nrow = 2))
      rownames(tab) <- c("0", "1")
      colnames(tab) <- c("0", "1")
      xtab <- caret::confusionMatrix(tab, positive = "1")
    }
  } else{
    xtab <- caret::confusionMatrix(table(as.factor(predY), as.factor(Y)), positive = "1")
  }
  xtab
}

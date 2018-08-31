############################################################################################
########################         Random Forests algorithm           ########################
############################################################################################

# =========================================================================================
# Impurity functions
# =========================================================================================
#' Title
#'
#' @param p 
#'
#' @return
#' @export
#'
#' @examples
entropy <- function(p) {
  if (any(p == 1))
    return(0)# works for the case when y has only 0 and 1 categorie...-sum(p *
  -sum(p *log(p, 2))
}

#' Title
#'
#' @param p 
#'
#' @return
#' @export
#'
#' @examples
gini <- function(p) {
  sum(p * (1 - p))
}
# =========================================================================================
# bsamples()
# Function which takes as inputs two integers "B" and "sampsize", a sample of data "data"
# and a boolean "replace".
# The function draws randomly B bootstrap samples of size "sampsize" (the size of the sample of data)
# and retunrs a matrix which with sampsize lines and B columns which contains the indices of
# the observationsbelonging each boostrap sample.
# =========================================================================================

#' Title
#'
#' @param ntree 
#' @param data 
#' @param sampsize 
#' @param replace 
#'
#' @return
#' @export
#'
#' @examples
bsamples <- function(ntree, data, sampsize, replace) {
  N <- nrow(data)
  bsamples <-
    apply(matrix(rep(1:N, ntree), ncol = ntree, nrow = N),
          2,
          sample,
          size = sampsize,
          replace = replace)
  oobsamples <- lapply(1:ntree, function(x)
    setdiff(1:N, bsamples[, x]))
  return(list(bsamples = bsamples, oobsamples = oobsamples))
}


# =========================================================================================
# group.selection()
# F
# =========================================================================================

#' Title
#' Function which takes inputs the vector "group" which indicates the group that contain each
#' variable and a number "mtry" (m<=length(group)).
#' The function selects randomly mtry groups among the length(group) groups and returns a
#' vector with the indices of the mtry selected group.
#' note that the default value for m is sqrt(length(group))
#' @param group 
#' @param mtry 
#'
#' @return
#' @export
#'
#' @examples
group.selection <-
  function(group, mtry = sqrt(unique(group[!is.na(group)]))) {
    if (length(which(is.na(group))) > 0) {
      print("Warning: there are NA in group; The fucntion will ignore these NA")
      label <- unique(group[!is.na(group)])
    } else{
      label <- unique(group)
    }
    p <- length(label)
    groups <- NULL
    if (mtry <= p) {
      groups <- sample(label, mtry, replace = FALSE)
    } else{
      print("WARNING: mtry must be lower or equal to p")
      groups <- label
    }
    return(groups)
  }



# =========================================================================================
# predict.rfgv()
# Functions which takes as inputs:
#     newdata         = a data frame containing the response value and the predictors and used to grow the tree
#     rfgvobject      = a vector with the group number of each variable
#
# =========================================================================================

# predict.rfgv<-function(newdata,rfgvobject){
#   pred<-matrix(rep(NA,nrow(newdata)*rfgvobject$ntree),nrow=nrow(newdata))
#   if(keep_forest==T){
#     if(rfgvobject$sampvar==TRUE & rfgvobject$sampvar_type=="2" & rfgvobject$maxdepth>1){
#       pred<-sapply(rfgvobject$forest,)
#     }else{
#
#     }
#   }else{
#     print("Prediction can be computed, restart the rfgv() function with the parameter keep_forest=TRUE")
#   }
#
# }

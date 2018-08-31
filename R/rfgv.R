#' Title Random Forest Grouped Variables
#'
#' @param data data frame containing the response value and the predictors and used to grow the tree
#' @param group a vector with the group number of each variable
#' (WARNING : if there are p goups, the groups must be numbers from 1 to p in increasing order)
#' @param groupImp a vector with the group number of each variable ==> the group used to compute the variable importance
#' @param ntree the number of trees to grow
#' @param mtry_group the number of variables randomly samples as candidates at each split
#' @param maxdepth only used with method="CARTGV". Integer indicating the maximal depth for a split-tree.
#'                   The default value is 2.
#' @param kfold only used with method="TPLDA". Integer indicating the numberof folds required
#'                   in the crossvalidation used to select the value of lambda.
#' @param replace Drawn of the bootstrap samples with or without replacment
#' @param sampsize the size of each boostrap sample
#' @param case_min minimun number of cases/non cases in a terminal nodes. For CARTGV, the default
#'                    is 1 while it is the number of fold for TPLDA.
#' @param grp.importance a boolean indicating the claculation or not of the importance of eacg group
#' @param test an independent data frame containing the same value that data
#' @param keep_forest a boolean indicating if the forest will be retained in the output object
#' @param crit the impurity function used (1=Gini/2=Entropie/3=Misclas)
#' @param penalty 
#' @param sampvar a boolean indicating if within each tree-split, a subset of variables 
#'                is drawn for each group
#' @param sampvar_type if sampvar_type="1" for each group a subset of variables belonging 
#'                     to the group is drawn at the beging of the
#'                     tree-split and only these variables are used to build the tree-split.
#'                     if sampvar_type="2" for each group, a subset of variables belonging to 
#'                     the group is drawn before each split withion the tree-split.
#' @param mtry_var usefull only if sampvar=TRUE. It indicates the number of drawn variables
#'
#' @return
#' @export
#'
#' @examples
rfgv <-
  function(data,
           group,
           groupImp,
           ntree = 200,
           mtry_group = floor(sqrt(length(unique(group[!is.na(group)])))),
           maxdepth = 1,
           kfold = 3,
           replace = T,
           sampsize = ifelse(replace == T, nrow(data), floor(0.632 * nrow(data))),
           case_min = 1,
           grp.importance = TRUE,
           test = NULL,
           keep_forest = F,
           crit = 1,
           penalty = "No",
           sampvar = FALSE,
           sampvar_type = "1",
           mtry_var = sapply(as.numeric(table(group[!is.na(group)])), function(x)
             floor(sqrt(x)))) {
    
    nbgrp <- length(unique(group[!is.na(group)]))
    predicted <- matrix(rep(NA, nrow(data) * ntree), ncol = ntree)
    vote <- matrix(rep(NA, nrow(data) * 2), ncol = 2)
    err <- rep(0, ntree)
    if (!is.null(test)) {
      vote.test <- matrix(rep(NA, nrow(test) * 2), ncol = 2)
      test_predicted <-
        matrix(rep(NA, nrow(test) * ntree), ncol = ntree)
      err.test <- rep(0, ntree)
    }
    if (keep_forest == TRUE) {
      forest <- list()
    } else{
      forest <- NULL
    }
    if (grp.importance == TRUE) {
      nbgrpImp <- length(unique(groupImp[!is.na(groupImp)]))
      acc <- rep(NA, 1)
      decreaacc <-
        matrix(rep(0, ntree * length(unique(groupImp[!is.na(groupImp)]))), nrow =
                 ntree)
      MeanDecrAcc <- rep(0, nbgrpImp)
    }
    samples <- bsamples(ntree = ntree, data, sampsize, replace)
    for (b in 1:ntree) {
      #if(sampvar==TRUE & sampvar_type=="2" & maxdepth>1){
      tree <-
        cartgv.rf(
          data[unlist(samples$bsamples[, b]),],
          group,
          crit = crit,
          case_min = case_min,
          maxdepth = maxdepth,
          p = mtry_group,
          penalty = penalty,
          mtry_var = mtry_var
        )
      predicted[unlist(samples$oobsamples[[b]]), b] <-
        as.numeric(as.character(
          predict.cartgv.rf(data[unlist(samples$oobsamples[[b]]),], tree$tree, tree$tree_split)$hat.Y
        ))
      
      #}else{
      #tree<-cartgv(data[unlist(samples$bsamples[,b]),],group,crit=crit,case_min=case_min,maxdepth=maxdepth,
      #            p=mtry_group,RF=T, IMPORTANCE=FALSE,penalty=penalty,sampvar=sampvar,mtry_var=mtry_var)
      #predicted[unlist(samples$oobsamples[[b]]),b]<-as.numeric(as.character(predict.cartgv(data[unlist(samples$oobsamples[[b]]),], tree$tree,tree$carts,tree$tables_coupures)$hat.Y))
      
      #}
      err[b] <-
        length(which(predicted[unlist(samples$oobsamples[[b]]), b] != data$Y[unlist(samples$oobsamples[[b]])])) /
        length(data$Y[unlist(samples$oobsamples[[b]])])
      if (grp.importance == TRUE) {
        if (sampvar == TRUE & sampvar_type == "2" & maxdepth > 1) {
          acc <-
            1 - as.numeric(impurity.cartgv.rf(data[unlist(samples$oobsamples[[b]]),], list(tree$tree), tree)$impurete$Misclass[1])
        } else{
          acc <-
            1 - as.numeric(impurity.cartgv(data[unlist(samples$oobsamples[[b]]),], list(tree$tree), tree)$impurete$Misclass[1])
        }
        groupselec <-
          unique(as.numeric(as.character(tree$tree$var[which(tree$tree$action == "1")])))
        if (sampvar == TRUE & sampvar_type == "2" & maxdepth > 1) {
          groupselecImp <- unique(groupImp[which(group %in% groupselec)])
          decreaacc[b, groupselecImp] <-
            as.numeric(unlist(sapply(groupselecImp,
                                     function(x)
                                       grpimpperm.rf(
                                         num_group = x,
                                         data = data,
                                         oobsamples = samples$oobsamples[[b]],
                                         group =
                                           groupImp,
                                         tree = tree,
                                         impurityacc = acc
                                       ))))
        } else{
          groupselecImp <- unique(groupImp[which(group %in% groupselec)])
          decreaacc[b, groupselecImp] <-
            as.numeric(unlist(sapply(groupselecImp,
                                     function(x)
                                       grpimpperm(
                                         num_group = x,
                                         data = data,
                                         oobsamples = samples$oobsamples[[b]],
                                         group =
                                           groupImp,
                                         tree = tree,
                                         impurityacc = acc
                                       ))))
        }
        acc <- NULL
      }
      if (!is.null(test)) {
        if (sampvar == TRUE & sampvar_type == "2" & maxdepth > 1) {
          test_predicted[, b] <-
            as.numeric(as.character(
              predict.cartgv.rf(test, tree$tree, tree$tree_split)$hat.Y
            ))
        } else{
          test_predicted[, b] <-
            as.numeric(as.character(
              predict.cartgv(test, tree$tree, tree$carts, tree$tables_coupures)$hat.Y
            ))
        }
        err[b] <-
          length(which(test_predicted[, b] != test$Y)) / length(test$Y)
        
      }
      if (keep_forest == TRUE) {
        #permet de faire ensuite de la prédiction
        forest[[b]] <- tree
      }
      rm(tree)
      gc()
    }
    cum.err <- apply(cbind(cumsum(err), 1:ntree), 1, function(x)
      x[1] / x[2])
    err.rate <-
      rowMeans(abs(predicted - as.numeric(as.character(data$Y))), na.rm = T)
    vote[, 2] <- rowMeans(predicted, na.rm = T)
    vote[, 2] <- ifelse(!is.nan(vote[, 2]), vote[, 2], NA)
    vote[, 1] <- ifelse(!is.na(vote[, 2]), 1 - vote[, 2], NA)
    vote <- as.data.frame(vote)
    names(vote) <- c("hat.Y=0", "hat.Y=1")
    pred <-
      ifelse(vote[, 1] > vote[, 2], "0", ifelse(vote[, 1] <= vote[, 2], "1", NA))
    oob.times <- rowSums(ifelse(!is.na(predicted), 1, 0), na.rm = T)
    confusion <-
      xtab_function(pred[!is.na(pred)], data$Y[!is.na(pred)])
    
    if (!is.null(test)) {
      cum.err.test <-
        apply(cbind(cumsum(err.test), 1:ntree), 1, function(x)
          x[1] / x[2])
      err.rate.test <-
        rowMeans(abs(test_predicted - as.numeric(as.character(test$Y))), na.rm =
                   T)
      vote.test[, 2] <- rowMeans(test_predicted, na.rm = T)
      vote.test[, 1] <-
        ifelse(!is.na(vote.test[, 2]), 1 - vote.test[, 2], NA)
      vote.test <- as.data.frame(vote.test)
      names(vote.test) <- c("hat.Y=0", "hat.Y=1")
      pred.test <-
        ifelse(vote.test[, 1] > vote.test[, 2],
               "0",
               ifelse(vote.test[, 1] <= vote.test[, 2], "1", NA))
      confusion.test <- xtab_function(pred.test, test$Y)
    } else{
      cum.err.test <- NULL
      err.rate.test <- NULL
      pred.test <- NULL
      confusion.test <- NULL
      vote.test <- NULL
    }
    if (grp.importance == TRUE) {
      MeanDecrAcc <- colSums(decreaacc) / ntree
      MeanDecrAccNor <-
        apply(cbind(MeanDecrAcc, as.numeric(table(groupImp[-1]))), 1, function(x)
          x[1] / x[2])
      importance <-
        as.data.frame(cbind(MeanDecrAcc, MeanDecrAccNor))
    } else{
      importance <- NULL
    }
    
    return(
      list(
        forest = forest,
        predicted = predicted,
        importance = importance,
        err.rate = err.rate,
        vote = vote,
        pred = pred,
        confusion = confusion,
        err.rate.test = err.rate.test,
        vote.test = vote.test,
        pred.test = pred.test,
        confusion.test = confusion.test,
        oob.times = oob.times,
        cum.err = cum.err,
        cum.err.test = cum.err.test,
        keep_forest = keep_forest,
        sampvar = sampvar,
        sampvar_type = sampvar_type,
        maxdepth = maxdepth,
        mtry_group = mtry_group,
        mtry_var = mtry_var,
        ntree = ntree
      )
    )
  }

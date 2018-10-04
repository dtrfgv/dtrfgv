############################################################################################
########################         COST COMPLEXITY PRUNING      ########################
############################################################################################

#### Followed the algorithm introduced by Breiman(page 293)


#' descendant
#'
#' function that find the direct descendants of a node
#'
#' @param node the number of the node in the tree \code{tree}.
#' @param tree a tree. A data frame returned by the function \code{CARTGV} and which summarizes the resulted CARTGV tree.
#'
#' @return descendants: a vector with the number of the descendants of \code{node} in \code{tree}.
#' @export
#'
descendant<-function(node,tree){
  descendants<-c()
  if(tree$leave[tree$node==node]=="*"){# if node is terminal
    descendants<-c(0,0)
  }else{
    descendants<-as.numeric(as.character(tree$node[which(tree$parent==node)]))
  }
  return(descendants)
}


#' conditional_error
#'
#' function that calculate the missclassification error assuming that we are in the node
#'
#' @param node the number of the node in the tree \code{tree}.
#' @param tree a tree. A data frame returned by the function \code{CARTGV} and which summarizes the resulted CARTGV tree.
#'
#' @return cond_err: an numeric indicating the missclassification error assuming that we are in the node.
#' @export
#'
conditional_error<-function(node,tree){
  cond_err<-NULL
  if(as.numeric(as.character(tree$yval[tree$node==node]))==1){
    cond_err<-(1-as.numeric(as.character(tree$prob[tree$node==node])))
  }else{
    cond_err<-as.numeric(as.character(tree$prob[tree$node==node]))
  }
  return(cond_err)
}


#' node_prob
#'
#' function that calculate the probability of being at the node
#'
#' @param node the number of the node in the tree \code{tree}.
#' @param tree a tree. A data frame returned by the function \code{CARTGV} and which summarizes the resulted CARTGV tree.
#'
#' @return a numeric indicating the probability that an observation falls into a node.
#' @export
#'
node_prob<-function(node,tree){
  return((as.numeric(as.character(tree$n[tree$node==node])))/(as.numeric(as.character(tree$n[1]))))
}



#' local_error
#'
#' function that calculate the missclassification error that the node contribute
#'
#' @param node the number of the node in the tree \code{tree}.
#' @param tree a tree. A data frame returned by the function \code{CARTGV} and which summarizes the resulted CARTGV tree.
#'
#' @return a numeric indicating the missclassification error that the node contribute.
#' @export
#'
local_error<-function(node,tree){
  err<-NULL
  if(as.numeric(as.character(tree$yval[tree$node==node]))==1){
    err<-as.numeric(as.character(tree$n_noncase[as.numeric(as.character(tree$node))==node]))/(as.numeric(as.character(tree$n[1])))
  }else{
    err<-as.numeric(as.character(tree$n_case[as.numeric(as.character(tree$node))==node]))/(as.numeric(as.character(tree$n[1])))
  }
  return(err)
}



#' get_nodes
#'
#'
#'function that gives the nodes of the subtree stremming from node
#'
#' @param node the number of the node in the tree \code{tree}.
#' @param tree a tree. A data frame returned by the function \code{CARTGV} and which summarizes the resulted CARTGV tree.
#'
#' @return a vector indicating the nodes of the subtree stremming from \code{node}.
#' @export
#'
get_nodes<-function(node,tree){
  nodes<-c(node)
  if(sum(descendant(node,tree))>0){
    nbChilds<-length(descendant(node,tree))
    for(child in 1:nbChilds){
      nodes<-c(nodes,get_nodes(descendant(node,tree)[child],tree))
    }
  }
  return(as.numeric(as.character(nodes)))
}


#' get_node_leaves
#'
#'function that gives the leaves of the subtree stremming from node
#'
#' @param node the number of the node in the tree \code{tree}
#' @param tree a tree. A data frame returned by the function \code{CARTGV} and which summarizes the resulted CARTGV tree.
#'
#' @return a vector indicating the leaves of the subtree stremming from \code{node}.
#' @export
#'
get_node_leaves<-function(node,tree){
  all_nodes<-get_nodes(node,tree)
  return(as.numeric(as.character(all_nodes[tree$leave[all_nodes]=="*"])))
}


#' node_error
#'
#' Function that calculate the misclassification error of the subtree stremming from node
#  Note that the error is not considered conditionaly to the fact that the root of the subtree is node
#'
#' @param node the number of the node in the tree \code{tree}.
#' @param tree a tree. A data frame returned by the function \code{CARTGV} and which summarizes the resulted CARTGV tree.
#'
#' @return a numeric indicating the misclassification error of the subtree stremming from node.
#' @export
#'
node_error<-function(node,tree){
  all_leaves<-get_node_leaves(node,tree)
  err<-0
  for(leave in all_leaves){
    err<-err+local_error(leave,tree)
  }
  return(err)
}


#' cost_error
#'
#' function that calculate the misclassification error of the subtree stremming from \code{node} penalized by alpha
#'
#' @param node the number of the node in the tree \code{tree}
#' @param tree a tree. A data frame returned by the function \code{CARTGV} and which summarizes the resulted CARTGV tree.
#' @param alpha a positive real. It is the value of the complexity paramter in the cost-complexity criterion.
#'
#' @return a numeric indicating the misclassification error of the subtree stremming from node penalized by alpha.
#' @export
#'
cost_error<-function(node,tree,alpha){
  return(node_error(node,tree)+ alpha*length(get_node_leaves(node,tree)))
}



#' node_cost
#'
#' function that gives the change in tree cost (difference between the node error and the local node error)  ==> gives alpha (function g(t) in the Breiman's book)
#'
#' @param node the number of the node in the tree \code{tree}
#' @param tree a tree. A data frame returned by the function \code{CARTGV} and which summarizes the resulted CARTGV tree.
#'
#' @return a numeric indicating the change in tree cost (difference between the node error and the local node error). It corresponds to the value of the complexity parameter \code{alpha}.
#' @export
#'
node_cost<-function(node,tree){
  return((local_error(node,tree)-node_error(node,tree))/(length(get_node_leaves(node,tree))-1))
}


#'best_node_cost
#'
#'
#' function that gives the minimum change in tree cost ==> gives alpha = min g(t), for all t in T
#'
#' @param node the number of the node in the tree \code{tree}
#' @param tree a tree. A data frame returned by the function \code{CARTGV} and which summarizes the resulted CARTGV tree.
#'
#' @return a numeric that gives the minimum change in the subtree stremming from \code{node}
#' @export
#'
best_node_cost<-function(node,tree){
  best_cost<-node_cost(node,tree)
  if(sum(descendant(node,tree))>0){# the node is not a leave
    nbChilds<-length(descendant(node,tree))
    for(child in 1:nbChilds){
      best_cost<-min(best_cost,best_node_cost(descendant(node,tree)[child],tree))
    }
  }else{# the node is a leave
    best_cost<-Inf
  }
  return(as.numeric(as.character(best_cost)))
}



#' cost_complexity
#'
#' @param tree a tree. A data frame returned by the function \code{CARTGV} and which summarizes the resulted CARTGV tree.
#'
#' @return a list containing the elements:
#'         cp: a data frame that summarizes the complexity as explained in the example in Breiman (chap 10, 1984). (k=number of iterations,
#'             (t_prunek=the k-th pruned node. All the subtree stremming from this node is deleted, N_Ttk= the number of leaves in the k-th
#'              pruned tree, R_Ttk=error for the pruned tree at the k-th step).
#'         summary: data frame indicating information about nodes in the maximal tree.
#'
#' @export
#'
cost_complexity<-function(tree){
  nb_noeuds<-length(tree$node)
  table_complexity<-as.data.frame(matrix(rep(NA,nb_noeuds*7),ncol=7,nrow=nb_noeuds))
  names(table_complexity)<-c("t","p_t","R_t","S_t","N_t","g_t","G_t")
  table_complexity$t<-as.numeric(as.character(tree$node))
  for(t in seq(nb_noeuds,1,-1)){
    if(sum(as.numeric(as.character(descendant(t,tree))))==0){
      table_complexity$p_t[t]<-as.numeric(as.character(tree$parent[as.numeric(as.character(tree$node[t]))]))
      table_complexity$N_t[t]<-1
      table_complexity$R_t[t]<-local_error(t,tree)
      table_complexity$S_t[t]<-node_error(t,tree)
      table_complexity$g_t[t]<-Inf
      table_complexity$G_t[t]<-Inf
    }else{
      table_complexity$p_t[t]<-as.numeric(as.character(tree$parent[as.numeric(as.character(tree$node[t]))]))
      table_complexity$N_t[t]<-length(get_node_leaves(t,tree))
      table_complexity$R_t[t]<-local_error(t,tree)
      table_complexity$S_t[t]<-node_error(t,tree)
      #table_complexity$g_t[t]<-node_cost(t,tree)
      table_complexity$g_t[t]<-(table_complexity$R_t[t]-table_complexity$S_t[t])/(table_complexity$N_t[t]-1)
      table_complexity$G_t[t]<-best_node_cost(t,tree)
    }
  }

  table<-table_complexity
  alpha<-0
  k<-c()
  N_Ttk<-c()
  alpha_k<-c()
  R_Ttk<-c()
  t_prunek<-c()
  e<-0.00005# conseill? par Breiman dans son livre
  t_prune<-NA # mettre NA car le premier sous arbre peut etre l'arbre maximal et dans ce cas aucun noeuds n'est retire
  g_t<-rep(NA,nb_noeuds)
  i<-1
  tab_g_k<-cbind(table_complexity$g_t)
  while(table$N_t[1]>1){
    ## Step 3
    if(table$G_t[1] > (alpha + e)){
      k<-c(k,i)
      N_Ttk<-c(N_Ttk,table$N_t[1])
      alpha_k<-c(alpha_k,alpha)
      R_Ttk<-c(R_Ttk,table$S_t[1])
      alpha<-table$G_t[1]
      t_prunek<-c(t_prunek,t_prune)
      t_prune<-NULL
      i<-i+1
    }

    if(table$N_t[1]>1){
      t=1
      while(table$G_t[t]<(table$g_t[t]-e)){
        childs=as.numeric(as.character(table$t[which(table$p_t==t)]))
        t<-childs[which(round(table$G_t[childs],6)==round(table$G_t[t],6))[1]]
      }

     if(is.null(t_prune)==TRUE || is.na(t_prune)){
        t_prune<-t
       }else{
         t_prune<-paste(t_prune,t,collapse="",sep=",")
       }


      ## Step 4
      table$N_t[t]<-1
      table$S_t[t]<- table$R_t[t]
      #table$g_t[t]<-Inf
      table$G_t[t]<-Inf
      ## Step 5
      while(t>1){
        t<-table$p_t[t]
        childs=as.numeric(as.character(table$t[which(table$p_t==t)]))
        table$N_t[t]<-table$N_t[childs[1]]
        table$S_t[t]<-table$S_t[childs[1]]
        G_t<-table$G_t[childs[1]]
        for(ichild in childs[-1]){
          table$N_t[t]<-table$N_t[t]+table$N_t[ichild]
          table$S_t[t]<-table$S_t[t]+table$S_t[ichild]
          G_t<-c(G_t,table$G_t[ichild])
        }
        table$g_t[t]<-(table$R_t[t]-table$S_t[t])/(table$N_t[t]-1)
        table$G_t[t]<-min(c(table$g_t[t],G_t))
      }
      tab_g_k<-cbind(tab_g_k,table$g_t)
    }
    # print(paste("k=",k))
    # print(paste("N_Ttk=",N_Ttk))
    # print(paste("alpha_k=",round(alpha_k,4)))
    # print(paste("R_Ttk=",round(R_Ttk,4)))
    # print(paste("t_prunek",t_prunek))
    # print(paste("alpha_courant",round(alpha,4)))
    # print(paste("t_prune_courant",t_prune))
    # print("********************************")
    # print("********************************")
    # print("********************************")
  }
  k<-c(k,i)
  N_Ttk<-c(N_Ttk,table$N_t[1])
  alpha_k<-c(alpha_k,alpha)
  R_Ttk<-c(R_Ttk,table$S_t[1])
  alpha<-table$G_t[1]
  t_prunek<-c(t_prunek,1)


  cp<-as.data.frame(cbind(k,N_Ttk,alpha_k=round(alpha_k,4),R_Ttk=round(R_Ttk,4),t_prunek))
  summary<-round(as.data.frame(cbind(table_complexity[,c(1,2,3,5)])),4)



  return(list(cp=cp,table_complexity=round(table_complexity,4),summary=summary,table_g_k=tab_g_k))
}



#' extratct_subtrees
#'
#' function that extract the sequence of nested trees
#'
#' @param tree a tree. A data frame returned by the function \code{CARTGV} and which summarizes the resulted CARTGV tree.
#' @param cp  a data frame that summarizes the complexity as explained in the example in Breiman (chap 10, 1984). It is the object \code{cp} returned by the function \code{cost_complexity}.
#'
#' @return a list containing the sequence of optimal subtrees.
#' @export
#'
extract_subtrees<-function(tree,cp){
  n<-nrow(cp)
  subtrees_seq<-list()
  tree_courant<-tree
  pruned_nodes_tot<-NULL
  for (i in 1:n) {
    if(stringr::str_detect(cp$t_prunek[i],",")!=TRUE || is.na(cp$t_prunek[i])){## si plusieurs sous-arbres sont retires a cette étape
      node<-cp$t_prunek[i]
    }else{
      node<-unlist(stringr::str_split(cp$t_prunek[i],","))
    }
    prune_nodes<-NULL
    for(j in 1:length(node)){
      if(!is.na(node[j])){
        prune_nodes<-c(prune_nodes,get_nodes(as.numeric(as.character(node[j])),tree_courant)[-1])
      }
    }

    if(length(prune_nodes)>0){
      subtree<-tree_courant[is.element(tree_courant$node,setdiff(tree_courant$node,prune_nodes)),]
      subtree$leave[which(as.numeric(as.character(subtree$node))%in%as.numeric(as.character(node)))]<-"*"
      subtree$action[which(as.numeric(as.character(subtree$node))%in%as.numeric(as.character(node)))]<--1
      subtree$var[which(as.numeric(as.character(subtree$node))%in%as.numeric(as.character(node)))]<-NA
      subtrees_seq[[i]] <- subtree
      tree_courant<-subtree
      subtree<-NULL
    }else{
      subtrees_seq[[i]] <- tree_courant
    }
  }
  return(subtrees_seq)
}


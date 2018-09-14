############################################################################################
########################         COST COMPLEXITY PRUNING      ########################
############################################################################################

#### Fait suivant les instructions de Breiman (page 293)

library(stringr)

# =========================================================================================
# descendant(node,tree) :  function that find the direct descendants of a node
# =========================================================================================

descendant<-function(node,tree){
  descendants<-c()
  if(tree$leave[tree$node==node]=="*"){# if node is terminal
    descendants<-c(0,0)
  }else{
    descendants<-as.numeric(as.character(tree$node[which(tree$parent==node)]))
  }
  return(descendants)
}

# =========================================================================================
# conditionnal_error(node,tree) : function that calculate the missclassification error 
# assuming that we are in the node 
# =========================================================================================

conditional_error<-function(node,tree){
  cond_err<-NULL
  if(as.numeric(as.character(tree$yval[tree$node==node]))==1){
    cond_err<-(1-as.numeric(as.character(tree$prob[tree$node==node])))
  }else{
    cond_err<-as.numeric(as.character(tree$prob[tree$node==node]))
  }
  return(cond_err)
}

# =========================================================================================
# node_prob(node,tree) : function that calculate the probability of being at the node
# =========================================================================================

node_prob<-function(node,tree){
  return((as.numeric(as.character(tree$n[tree$node==node])))/(as.numeric(as.character(tree$n[1]))))
}

# =========================================================================================
# local_error(node,tree) : function that calculate the missclassification error 
# that the node contribute
# WARNING : don't use the function local_error_0(node,tree) which use to much rounded number 
# and so lead to error
# =========================================================================================

local_error_0<-function(node,tree){
  return(node_prob(node,tree)*conditional_error(node,tree))
}


local_error<-function(node,tree){
  err<-NULL
  if(as.numeric(as.character(tree$yval[tree$node==node]))==1){
    err<-as.numeric(as.character(tree$n_noncase[as.numeric(as.character(tree$node))==node]))/(as.numeric(as.character(tree$n[1])))
  }else{
    err<-as.numeric(as.character(tree$n_case[as.numeric(as.character(tree$node))==node]))/(as.numeric(as.character(tree$n[1])))
  }
  return(err)
}



# =========================================================================================
# get_nodes(node,tree) : function that gives the nodes of the subtree stremming 
# from node
# WARNINGS: this function was modified on the 4th december 2017 to be adapted to 
#           non-binary trees.
# =========================================================================================

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

# =========================================================================================
# get_node_leaves(node,tree) : function that gives the leaves of the subtree stremming 
# from node
# =========================================================================================

get_node_leaves<-function(node,tree){
  all_nodes<-get_nodes(node,tree)
  return(as.numeric(as.character(all_nodes[tree$leave[all_nodes]=="*"])))
}

# =========================================================================================
# node_error(node,tree) : function that calculate the misclassification error of the subtree 
# stremming from node
# Note that the error is not considered conditionaly to the fact that the root of 
# the subtree is node
# =========================================================================================

node_error<-function(node,tree){
  #subtree<-as.data.frame(tree[which(tree$node %in% get_nodes(node,tree)),])
  all_leaves<-get_node_leaves(node,tree)
  err<-0
  for(leave in all_leaves){
    err<-err+local_error(leave,tree)
  }
  return(err)
}

# =========================================================================================
# cost_error(node,tree) : function that calculate the misclassification error of the subtree 
# stremming from node penalized by alpha
# =========================================================================================

cost_error<-function(node,tree,alpha){
  return(node_error(node,tree)+ alpha*length(get_node_leaves(node,tree)))
}


# =========================================================================================
# node_cost (node,tree) : function that gives the change in tree cost (difference between
# the node error and the local node error) 
# ==> gives alpha (function g(t) in the Breiman's book)
# =========================================================================================

node_cost<-function(node,tree){
  return((local_error(node,tree)-node_error(node,tree))/(length(get_node_leaves(node,tree))-1))
}

# =========================================================================================
# best_node_cost (node,tree) : function that gives the minimum change in tree cost 
# ==> gives alpha = min g(t), for all t in T
# WARNINGS: this function was modified on the 4th december 2017 to be adapted to 
#           non-binary trees.
# =========================================================================================

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


# =========================================================================================
# Cost_complexity (tree) :  function that built the sequence of alpha and of nested trees
# 
# =========================================================================================


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
## La fonction retourne:
## CP: tableau de complexite comme dans le livre de Breiman
##### k : numero de l'iteration
##### t_prunek : k-ieme noeuds elague 
##### (il devient une feuille, i.e. tous les neouds partant du sous-arbre issue de t_pruneK sont supprimés)
##### N_Ttk : nombre de feuilles du k-ieme arbre elague
##### alpha_k : valeur de la temperature a l'etape k
##### R_Ttk   : erreur de la partie elaguee a l'etape k

## table_complexity: tableau donnat des infos sur chaque noeuds de l'arbre maximal
##### k : numero de l'iteration
##### t_prunek : k-ieme noeuds elague (il devient une feuille)
##### N_Ttk : nombre de feuilles deu k-ieme arbre ?lagu?
##### alpha_k : valeur de la temperature a l'etape k
##### R_Ttk   : erreur de la partie elaguee a l'etape k



cost_complexity_binary<-function(tree){
  nb_noeuds<-length(tree$node)
  table_complexity<-as.data.frame(matrix(rep(NA,nb_noeuds*9),ncol=9,nrow=nb_noeuds))
  names(table_complexity)<-c("t","l_t","r_t","p_t","R_t","S_t","N_t","g_t","G_t")
  table_complexity$t<-as.numeric(as.character(tree$node))
  table_complexity$l_t<-rep(0,nb_noeuds)
  table_complexity$r_t<-rep(0,nb_noeuds)
  for(t in seq(nb_noeuds,1,-1)){
    if(length(as.numeric(as.character(descendant(t,tree))))==0){
      table_complexity$p_t[t]<-as.numeric(as.character(tree$parent[as.numeric(as.character(tree$node[t]))]))
      table_complexity$N_t[t]<-1
      table_complexity$R_t[t]<-local_error(t,tree)
      table_complexity$S_t[t]<-node_error(t,tree)
      table_complexity$g_t[t]<-Inf
      table_complexity$G_t[t]<-Inf
    }else{
      table_complexity$p_t[t]<-as.numeric(as.character(tree$parent[as.numeric(as.character(tree$node[t]))]))
      table_complexity[t,2:3]<-as.numeric(as.character(descendant(t,tree)))
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
        if(table$G_t[t]==table$G_t[table$l_t[t]]){
          t<-table$l_t[t]
        }else{
          t<-table$r_t[t]
        }
      }
      #g_t[t]<-alpha
      
      if(is.null(t_prune)==TRUE || is.na(t_prune)){
        t_prune<-t
      }else{
        t_prune<-paste(t_prune,t,collapse="",sep=",")
      }
      
      
      ## Step 4
      table$N_t[t]<-1
      table$l_t[t]<-0
      table$r_t[t]<-0
      table$S_t[t]<- table$R_t[t]
      #table$g_t[t]<-Inf
      table$G_t[t]<-Inf
      ## Step 5
      while(t>1){
        t<-table$p_t[t]
        table$N_t[t]<-table$N_t[table$l_t[t]]+table$N_t[table$r_t[t]]
        table$S_t[t]<-table$S_t[table$l_t[t]]+table$S_t[table$r_t[t]]
        table$g_t[t]<-(table$R_t[t]-table$S_t[t])/(table$N_t[t]-1)
        table$G_t[t]<-min(table$g_t[t],table$G_t[table$l_t[t]],table$G_t[table$r_t[t]])
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
## La fonction retourne:
## CP: tableau de complexite comme dans le livre de Breiman
##### k : numero de l'iteration
##### t_prunek : k-ieme noeuds elague 
##### (il devient une feuille, i.e. tous les neouds partant du sous-arbre issue de t_pruneK sont supprimés)
##### N_Ttk : nombre de feuilles du k-ieme arbre elague
##### alpha_k : valeur de la temperature a l'etape k
##### R_Ttk   : erreur de la partie elaguee a l'etape k

## table_complexity: tableau donnat des infos sur chaque noeuds de l'arbre maximal
##### k : numero de l'iteration
##### t_prunek : k-ieme noeuds elague (il devient une feuille)
##### N_Ttk : nombre de feuilles deu k-ieme arbre ?lagu?
##### alpha_k : valeur de la temperature a l'etape k
##### R_Ttk   : erreur de la partie elaguee a l'etape k


# =========================================================================================
# extratct_subtrees (tree,cp) :  function that extract the sequence of nested trees
# 
# =========================================================================================

extract_subtrees<-function(tree,cp){
  n<-nrow(cp)
  subtrees_seq<-list()
  tree_courant<-tree
  pruned_nodes_tot<-NULL
  for (i in 1:n) {
    if(str_detect(cp$t_prunek[i],",")!=TRUE || is.na(cp$t_prunek[i])){## si plusieurs sous-arbres sont retires a cette étape
      node<-cp$t_prunek[i]
    }else{
      node<-unlist(str_split(cp$t_prunek[i],","))
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


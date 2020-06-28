Eliminate<-function(tre,phylotable,parent=0,recursion=FALSE){
  #Given a phylogenetic tree and a corresponding phylotable, remove uninformative branches 
  #of the "parent" node from the phylotable. Uninformative branches are child branches 
  #having the same representative taxonomies with the "parent" node.
  #a litlle differnt to the function definition 
  
  
  #tre: phylogentic tree, datatype: phylo
  #phylotable: dataframe with a column named "Node"(numeric) and a column named "represent", 
  #            datatype: data.frame|matrix
  #parent: the parent node, defualt 0 means eliminating from the root of the tree,
  #        datatype: interger
  
  #value: a data.frame similar to phylotable but uninformative branch removed 
  
  #data type check
  if(!recursion){#first call
    if(!is.data.frame(phylotable)) stop("please check object \"phylotable\" is of class data.frame")
    phylotable[order(phylotable$Node),]->phylotable
    
    if(!inherits(tre,"phylo")) stop("object \"tre\" is not of class \"phylo\"")
    N=length(tre$tip.label)
    M=N+tre$Nnode

    if(!isTRUE(all.equal(phylotable$Node,1:M)))
      stop("object \"phylotable$Node\" and the \"tre\" doesn't match, 
           it may result from duplicated values in \"phylotable$Node\",
           or tips/nodes in \"tre\" doesn't match \"phylotable$Node\"")

    if(parent<1){ #setting parent to root
      require(ape)
      if(!ape::is.rooted(tre)) stop("object \"tre\" is not rooted")
      if(tre$node.label[1]=="Root") parent=N+1 #how to find root of a tree?
      else if(length(which(tre$node.label=="Root"))!=0) 
        parent=N+which(tre$node.label=="Root")
      else stop("cannot find the \"Root\" of the \"tre\"")
    }
  }
  
  N=length(tre$tip.label)
  M=N+tre$Nnode
  uninformative<-c()
  
  if(phylotable$represent[parent]=="NA"){
    uninformative<-c(uninformative,parent) #if the parent is not available, delete it
    if(parent>N&parent<=M){  #if parent is a node to find the children
      tre$edge[which(tre$edge[,1]==parent),2]->children
      if(length(children)!=2) warning("node with branch(s) that are not equal to 2")
      for(i in 1:length(children)) 
        uninformative=c(uninformative,Eliminate(tre,phylotable,children[i],recursion=TRUE))
    }
    else if(parent<=N) {} #if parent is a leaf, just itself
    else stop("the \"parent\" is out of the \"tre\"")
  }
  else{
    if(parent>N&parent<=2*N-1){
      uninformative=getDescendants(tre,parent) #if the parent is a node but assained, delete the descendants
    }
    else if(parent<=N) {}  #if the node is a leaf and assained, reserve it
    else stop("the \"parent\" is out of the \"tre\"")
  }
  
  if(recursion)  return(uninformative) else{ #first call
      return(phylotable[-uninformative,])
    }
}
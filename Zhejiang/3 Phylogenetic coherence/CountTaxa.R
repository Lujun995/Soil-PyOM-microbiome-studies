CountTaxa<-function(tre,taxonomy,parent=0,recursion=FALSE){
  #To count taxanomy composition of parent node/tip
  #tre: phylogenetic tre, datatype="phylo"
  #taxonmy: taxonomy information of tree leavea, datatype="data.frame", the first row is tiplabels 
  #         and the second row is their taxonomies.
  #parent: the node, datatype=interger, if==0, from the root of the tre.
  #recursion: flag for recursion usage
  
  #Values:
  #A list contains:
  #TaxaCount: the taxanomy composition of parent node/tip
  #PhyloTable: the taxanomy composition of parent as well as child node/tip
  
  #check type errors in the first call
  if(!recursion){ #first call
    if (!inherits(tre, "phylo")) stop("object \"tre\" is not of class \"phylo\"")
    if(parent<1){ #setting parent to root
      require(ape)
      if(!ape::is.rooted(tre)) stop("object \"tre\" is not rooted")
      N=length(tre$tip.label)
      if(tre$node.label[1]=="Root") parent=N+1 #how to find root of a tree?
      else if(length(which(tre$node.label=="Root"))!=0) 
        parent=N+which(tre$node.label=="Root")
      else stop("cannot find the \"Root\" of the \"tre\"")
    }
    
    if(length(which(!(taxonomy[,1]%in%tre$tip.label)))>0)
      warning("\"taxonomy\" contains more tiplab(s)")
    if(length(which(!(tre$tip.label%in%taxonomy[,1])))>0) 
      #maybe just children of the parent node in the taxonomy is OK
      stop("\"tre\" contains unknown tiplab(s)") #remove unused taxonomy levels and tree tips in the future
    
    if(inherits(taxonomy[,2],"character")) as.factor(taxonomy[,2])->taxonomy[,2]
    if(inherits(taxonomy[,2],"factor")) {
      length(levels(taxonomy[,2]))->n
      levels(taxonomy[,2])->levels.taxonomy
      as.numeric(taxonomy[,2])->taxonomy[,2]
    } else stop("data type of \"taxonomy\" is unkown")
  } else 
    n<-recursion
  
  N=length(tre$tip.label)
  M=N+tre$Nnode
  taxacount<-rep(0,times=n)
  phylolist<-c()
  if(parent<=N) {
    temp<-taxonomy[which(taxonomy[,1]==tre$tip.label[parent]),2]
    taxacount[temp]=taxacount[temp]+1
    phylolist<-c(parent,taxacount)
  }
  else if(parent>N&parent<=M){
    tre$edge->edge
    edge[which(edge[,1]==parent),2]->children
    if(length(children)==2){
      CountTaxa(tre,taxonomy,children[1],n)->a
      CountTaxa(tre,taxonomy,children[2],n)->b
      taxacount=a[[1]]+b[[1]]
      phylolist<-c(c(parent,taxacount),a[[2]],b[[2]])
    }else {
      warning("node with branch(s) that are not equal to 2")
      for(i in 1:length(children)) {
        CountTaxa(tre,taxonomy,children[i],n)->a
        taxacount=taxacount+a[[1]]
        phylolist<-c(phylolist,a[[2]])
      }
      phylolist<-c(c(parent,taxacount),phylolist)
    }
  }else stop("the \"parent\" is out of the \"tre\"")
  
  
  if(recursion) return(list(taxacount,phylolist)) else
  { #first call
    phylotable<-matrix(data=phylolist,ncol=n+1,byrow=TRUE)
    colnames(phylotable)<-c("Node",levels.taxonomy)
    phylotable[order(phylotable[,1]),]->phylotable
    names(taxacount)<-levels.taxonomy
    if(sum(duplicated(phylotable[,1]))>0)
      warning("a duplicated node in the phylotable, please check object \"tre\"!")
    return(list(taxacount,phylotable))
  }
}
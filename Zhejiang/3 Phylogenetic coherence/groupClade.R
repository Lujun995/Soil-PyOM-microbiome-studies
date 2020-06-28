groupClade<-function(tre,node,represent=NULL){
  #To clade the the tre
  #tre: a phylogenetic tre, datatype=phylo
  #node: a node list, datatype=numeric|interger
  #represent: a informative vector containing corresponding taxonomies of each node
  
  #values: a phylogenetic tre with a attribute list "group" 
  
  if(!inherits(tre,"phylo")) stop("object \"tre\" is not of class \"phylo\"") 
  if(is.null(represent)) represent=1:length(node)
  if(length(represent)!=length(node))
    stop("length of object \"represent\" is not equal to \"node\"")
  group<-rep("NA",time=ape::Nnode(tre)+ape::Ntip(tre))
  as.character(unique(represent))->taxa
  for(i in 1:length(taxa)) {
    getDescendants(tre,node[represent==taxa[i]])->childnodes
    childnodes<-childnodes[!duplicated(childnodes)]
#    childnodes-ape::Ntip(tre)->childnodes
    group[c(node[represent==taxa[i]],childnodes)]<-taxa[i]
  }
  as.factor(group)->group
  attributes(tre)$group<-group
  return(tre)
}
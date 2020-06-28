AnnotateTre<-function(tre,taxonomy,...){
  #To find taxonomy clade in the tre
  #tre: phylogenetic tree, datatype=phylo
  #taxonomy: taxonomy information table with the first row is OTU labels 
  #         and the left row(s) are their taxonomies.
  #...: other parameters passed to function ????how to pass ...to sub functions???
  
  #Other parameters as: parent,TOLERANT,AMBIGUITY,label ann Missing
  
  #value: a phylogenetic tree with attribute "group" and the node hubs for each taxonomies
  
  require(ape)
  require(ggtree)
  
  #datatype check
  if (!inherits(tre, "phylo"))
    stop("object \"tre\" is not of class \"phylo\"")
  if(!(class(taxonomy)%in%c("data.frame","matrix")))
    stop("object \"taxonomy\" is of undefined class")
  if(!ape::is.rooted(tre)) {
    warning("The tree given is unrooted. Now root it.")
    ape::root(tre,1,r=TRUE)->tre
  }
  CountTaxa(tre = tre,taxonomy = taxonomy[,1:2])[[2]]->phylotable
  AnnotationNode(phylotable = phylotable)->phylotable
  Eliminate(tre,phylotable)->cladehubs
#  cladehubs[which(cladehubs$Node>length(tre$tip.label)),]->cladehubs
  groupClade(tre,cladehubs$Node,cladehubs$represent)->tretemp 

  return(list(tretemp,cladehubs))

}




# as.phylo.formula.new(x,data=parent.frame(),node.label=TRUE,...){
#   as.phylo.formula(x,data,...)->phylo
#   if(node.label)
#   {
#     for(j in 6:1){
#       phylo$node.label<-rep("NA",Nnode(phylo))
#       taxa<-unique(data[j+1])
#       if(length(temp)==0) temp<-rep(NA,length(taxa))
#       for(i in 1:length(taxa)){
#         phylo$edge[phylo$edge[,2]==i,]->nodenum->temp
#         data[data[,j+1]==phylo$tip.label[i],j]->nodelabel
#         node.label[nodenum]<-nodelabel
#       }
#     temp<-nodenum
#     }
#   }  else {
#     
#   }
# }
# 
# node.label1<-test$tip.label
# node.num1<-1:Ntip(test)
# edges<-test$edge
# rownames(edges)<-edges[,2]
# test$node.label<-rep("NA",Nnode(test))
# 
# rownames(anno.derep)<-anno.derep$Genus
# anno.derep<-anno.derep[node.label1,]
# node.label2<-anno.derep$Family
# node.num2<-edges[as.character(node.num1),1]
# test$node.label[node.num2]<-node.label2
# 
# node.label2[!duplicated(node.label2)]->node.label1
# node.num2[!duplicated(node.num2)]->node.num1
# anno.derep<-anno.derep[!duplicated(anno.derep$Family),]
# 
# rownames(anno.derep)<-anno.derep$Family
# anno.derep<-anno.derep[node.label1,]
# node.label2<-anno.derep$Order
# node.num2<-edges[as.character(node.num1),1]
# test$node.label[node.num2]<-node.label2
# 
# node.label2->node.label1
# node.num2->node.num1
# anno.derep<-anno.derep[!duplicated(anno.derep$Order),]
# 
# taxa=unique(anno.derep$Domain)
# for(i in 1:length(taxa)){
#   tip=anno.derep$Genus[anno.derep$Domain==taxa[i]]
#   cat("Node=",getMRCA(test,tip),"\n")
# }

as.phylo.data.frame<-function(data){
  ncol(data)->n
  Ntip<-length(unique(data[,n]))
  Nnode=1
  for(i in 1:(n-1)){
    Nnode=Nnode+length(unique(data[,i]))
  }

  seque<-paste0(paste0(rep("data[,",n),1:n,rep("]",n)),collapse = ",")
  eval(parse(text=paste("data<-data[order(",seque,"),]")))
  for(i in 1:n) data[,i]<-as.character(data[,i])

#  node.label<-rep("",Nnode)
  node.label<-"root"
  tip.label<-data[,n]
  edge<-NULL
  for(i in 1:(n-1)){
    taxa<-unique(data[,i])
    if(i==1){
      top<-1
      c(node.label,taxa)->node.label
      length(taxa)->k
      end<-top+length(taxa)
      edge<-matrix(c(rep(1,k),1+1:k),ncol=2)
    }
    for(j in 1:length(taxa)){
      childtaxa<-unique(data[data[,i]==taxa[j],i+1])
      length(node.label)->end
      length(childtaxa)->k
      c(node.label,childtaxa)->node.label
      edge<-rbind(edge,matrix(c(rep(top+j,k),end+1:k),ncol=2))
    }
    top<-top+j
  }
  edge=edge+Ntip
  edge[edge>Nnode+Ntip]<-edge[edge>Nnode+Ntip]-Ntip-Nnode
  node.label<-node.label[1:Nnode]
  phylo<-list(edge=edge,Nnode=Nnode,tip.label=tip.label,node.label=node.label)
  attributes(phylo)$class<-"phylo"
  attributes(phylo)$order<-"cladewise"
  return(phylo)
}



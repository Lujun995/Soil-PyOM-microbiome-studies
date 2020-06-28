library(ggtree)
library(ape)
library(picante)
library(caper)
source('comparative.data.R')
source('phylo.d.R')
source('summary.cladehub.R')
source("SigClade4.R")
source("as.phylo.data.frame.R")
source('pSpecialists.R')


dget("anno.data")->anno.data
anno.data$response->response
response[!anno.data$ifsignificant]<-0
response<-as.factor(response)
levels(response)<-c("Positive","NA","Negative")
anno.data$response<-response
anno.data<-anno.data[,-ncol(anno.data)]
response<-as.data.frame(response)
rownames(response)<-rownames(anno.data)
#Significant taxa (not neccessarily conservatism part!!)
#SigClade2(anno.data[,1:7],anno.data[anno.data$response=="Positive",1:7])->PyOMPhyloNew
# SigClade3(anno.data[,1:7],which(anno.data$response=="Positive"),Bonferroni=0.05)->PyOMPhyloNewNew
#SigClade2(anno.data[,1:7],anno.data[anno.data$response=="Negative",1:7])->SoilPhyloNew
#SigClade3(anno.data[,1:7],which(anno.data$response=="Negative"),Bonferroni=0.05)->SoilPhyloNewNew
SigClade4(anno.data[,1:7],which(anno.data$response=="Positive"),which(anno.data$response=="Negative"),Bonferroni=0.05)->all.coherent
temp<-matrix(0, nrow=nrow(all.coherent),ncol=3)
colnames(temp)<-c("PyOM","Neutral","Soil")
cbind(all.coherent,temp)->all.coherent
for(i in 1:nrow(all.coherent)) {
  anno.data[,c(as.character(all.coherent$Phylogeny)[i],"response")]->temp
  summary(temp[temp[,1]==as.character(all.coherent$taxa)[i],2])[c("PyOM","NA's","Soil")]->all.coherent[i,5:7]
}
write.csv(all.coherent,"CoherenceUnits.csv")

#Coherence test 
#not run
# P=rep(1,dim(SoilPhyloNewNew)[1])
# counts=matrix(rep(0,2*dim(SoilPhyloNewNew)[1]),ncol=2)
# for(i in 1:dim(SoilPhyloNewNew)[1]) try({
#   subtreeOTU=as.character(anno.data[anno.data[,as.character(SoilPhyloNewNew$Phylogeny)[i]]==as.character(SoilPhyloNewNew$taxa)[i],1])
#   dropping=setdiff(as.character(annotatedtre[[1]]$tip.label),subtreeOTU)
#   subTree=drop.tip(annotatedtre[[1]],dropping)
#   phylo.d(category,subTree,"OTUID","negatives",permut = 9999)->temp
#   temp$Pval1->P[i]
#   temp$StatesTable->counts[i,]
# })
# P->OverDispersal
# colnames(counts)<-c("NullNum","SigNum")
# cbind(SoilPhyloNewNew,OverDispersal,counts)->SoilPhyloNewNew
# 
# P=rep(1,dim(PyOMPhyloNewNew)[1])
# D=rep(1,dim(PyOMPhyloNewNew)[1])
# counts=matrix(rep(0,2*dim(PyOMPhyloNewNew)[1]),ncol=2)
# for(i in 1:dim(PyOMPhyloNewNew)[1]) try({
#   subtreeOTU=as.character(anno.data[anno.data[,as.character(PyOMPhyloNewNew$Phylogeny)[i]]==as.character(PyOMPhyloNewNew$taxa)[i],1])
#   dropping=setdiff(as.character(annotatedtre[[1]]$tip.label),subtreeOTU)
#   subTree=drop.tip(annotatedtre[[1]],dropping)
#   phylo.d(category,subTree,"OTUID","positives",permut = 9999)->temp
#   temp$Pval1->P[i]
#   temp$StatesTable->counts[i,]
# })
# P->OverDispersal
# colnames(counts)<-c("NullNum","SigNum")
# cbind(PyOMPhyloNewNew,OverDispersal,counts)->PyOMPhyloNewNew



# rbind(SoilPhyloNewNew,PyOMPhyloNewNew)->CoherenceUnits
# Type<-c(rep("Soil",dim(SoilPhyloNewNew)[1]),rep("PyOM",dim(PyOMPhyloNewNew)[1]))
# cbind(CoherenceUnits,Type)->CoherenceUnits
# CoherenceUnits[,"OverDispersal"]<-p.adjust(CoherenceUnits[,"OverDispersal"],method="bonferroni")
# write.csv(CoherenceUnits,"CoherenceUnits.csv")


#Visualization in the phylogenetic tree, see SnowFlake.R


#Visualization
# anno.drop<-anno.data[!(anno.data$Phylum=="Unassigned"),]
# anno.drop<-anno.drop[!(anno.data$Genus==""),]
# anno.derep<-anno.drop[!duplicated(anno.drop$Genus),]
# c(test2$tip.label,test2$node.label)->test
# as.phylo.data.frame(anno.derep[,2:7])->test2
# group=rep(0,length(test))
# for(i in 1:nrow(PyOMPhyloNewNew)){
#   which(test==PyOMPhyloNewNew$taxa[i])->node.num
#   getDescendants(test2,node.num)->childnodes
#   group[c(node.num,childnodes)]=group[c(node.num,childnodes)]+1
# }
# for(i in 1:nrow(SoilPhyloNewNew)){
#   which(test==SoilPhyloNewNew$taxa[i])->node.num
#   getDescendants(test2,node.num)->childnodes
#   group[c(node.num,childnodes)]=group[c(node.num,childnodes)]-1
# }
# 
# index=(test %in% c(as.character(PyOMPhyloNewNew$taxa),as.character(SoilPhyloNewNew$taxa)))
# index2=rep(0,length(index))
# index2[test %in% as.character(PyOMPhyloNewNew$taxa)]=
#   index2[test %in% as.character(PyOMPhyloNewNew$taxa)]+3
# index2[test %in% as.character(SoilPhyloNewNew$taxa)]=
#   index2[test %in% as.character(SoilPhyloNewNew$taxa)]-4
# 
# ggtree(test2,layout="fan",aes(color=group)) +
#   scale_colour_gradient2(low="firebrick3", mid="#D4D4D4",high="green4")+
#   geom_point2(aes(shape=isTip, color=index2), size=3)+
#   geom_point2(aes(subset=index,shape=isTip, color=index2), size=3)+
#   geom_text2(aes(subset=index,label=test,colour=index2),
#              position=position_jitter(width=0.7,height=0.7),fontface =2)+
#   theme(legend.position="right")








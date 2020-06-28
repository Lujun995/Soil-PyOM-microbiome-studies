library(ggtree)
source("gheatmap.R")
source("AnnotateTre.R")
source("groupClade.R")
source("CountTaxa.R")
source("Eliminate.R")
source("getDescendants.R")
source("AnnotationNode.R")
#also dependent on the result from SignificantTaxonomy.R

# ggtree(root.phylotre,layout="radial",color="#D4D4D4")+xlim(0, 1.8)->p
# p %<+% anno + 
#   geom_tippointtest(aes(size=size, color=corr), alpha=0.8)+
#   theme(legend.position="right")+scale_color_manual(values = c(color3,color1))->p1
# p1
# 
# 
# ggtree(temp[[1]],layout="fan",color="#D4D4D4")->p
# levels(temp[[2]]$represent)->representlevel
# length(representlevel)->ni
# for(i in 1:ni){
#   temp[[2]]$Node[which(temp[[2]]$represent==representlevel[i])]->subnode
#   for(j in 1:length(subnode))
#     p=p+geom_hilight(subnode[j], fill=rainbow(ni)[i])
# }
# p+xlim(0, 1.8)+theme(legend.position="right")->p
# p
# 
# AnnotateTre(refinedtre,anno.data.test[,c(1,6)])->temp
# ggtree(temp[[1]],layout="fan",aes(color=group))+theme(legend.position="right")->p
# p %<+% anno.data.test + geom_tippoint(na.rm = TRUE,aes(color=corr), alpha=0.8)
# 
# AnnotateTre(refinedtre,anno.data.test[,c(1,6)])->temp
# ggtree(temp[[1]],layout="fan",aes(color=group))+theme(legend.position="right")->p
# p %<+% anno + geom_tippoint(na.rm = TRUE,aes(color=corr), alpha=0.8)

load("D:\\World Forest Data\\Fuyang\\data\\.Rdata")
dget("anno.data")->anno.data
anno.data$response->response
response[!anno.data$ifsignificant]<-0
response<-as.factor(response)
levels(response)<-c("PyOM",NA,"Soil")  #change level -1 to "Positive" to PyOM addition
anno.data$response<-response
anno.data<-anno.data[,-ncol(anno.data)]
response<-as.data.frame(response)
rownames(response)<-rownames(anno.data)

all.coherent[order(all.coherent$Phylogeny),]->temp3
anno.data[,as.character(temp3$Phylogeny)]->temp2
coherence<-rep(NA,dim(anno.data)[1])
for(i in ncol(temp2):1){
  coherence[which(as.character(temp2[,i])==as.character(temp3$taxa[i]))]<-as.character(temp3$type[i])
}
rm(temp2,temp3)
cbind(as.factor(coherence),response)->response
factor(response$response,levels = c("PyOM","Soil"))->response$response
levels(response$`as.factor(coherence)`)[2]<-"Soil"

AnnotateTre(root.phylotre,anno.data[,c(1,3)])->annotatedtre
attributes(annotatedtre[[1]])$group->group
as.character(group)->group
unique(group)->taxa
for(i in 1:length(taxa)) if(length(group[group==taxa[i]])<61) group[group==taxa[i]]<-"Other"
as.factor(group)->attributes(annotatedtre[[1]])$group
ggtree(annotatedtre[[1]],layout="fan",aes(color=group))+theme_minimal()+
  theme(panel.grid = element_blank(),legend.position="right",
        axis.text.x = element_blank(),axis.text.y = element_blank())->p
#+xlim(0, 1.8)->p
p %<+% anno.data + #geom_tippoint(aes(color=corr,subset=test), alpha=1,size=2) +
  scale_color_manual(values = c(tol18rainbow[1:11],
                                "#D4D4D4",tol18rainbow[12:16],
                                "#888888",tol18rainbow[17:18]))->p
gheatmap(p,response,colnames = FALSE,width = 0.05,offset = -2)->p1
p1+scale_fill_manual(values=c("green4","firebrick3"),na.value = NA)
#rotate_tree(p, 90)
 
# AnnotateTre(refinedtre,anno.data.test[,c(1,10)])->temp
# ggtree(temp[[1]],layout="fan",aes(color=group))+theme(legend.position="right")->p
# p %<+% anno.data.test + geom_tippoint(aes(color=corr,subset=test), alpha=0.8)
# 
# AnnotateTre(refinedtre,anno.data.test[,c(1,6)])->temp
# ggtree(temp[[1]],layout="fan",aes(color=group))+theme(legend.position="right")->p
# p %<+% anno.data.test + geom_tippoint(aes(color=corr,subset=test), alpha=0.8,size=2) +
#   scale_color_manual(values = c(color3,color1,rainbow(39,s=0.6,v=0.75,alpha=0.8)[1:22],"#D4D4D4",
#                                 rainbow(39,s=0.6,v=0.75,alpha=0.8)[23:36],"#888888",
#                                 rainbow(39,s=0.6,v=0.75,alpha=0.8)[37:39]))->p
# temp<-list();for(i in 1:6) AnnotateTre(refinedtre,anno.data.test[,c(1,i+4)])->temp[[i]]
# phylo2node<-temp[[1]][[2]][,c(1,ncol(temp[[1]][[2]]))]
# for(i in 2:6) rbind(phylo2node,temp[[i]][[2]][,c(1,ncol(temp[[i]][[2]]))])->phylo2node
# phylo2node[phylo2node[,2] %in% SigSoilTaxa$Name,]->SigSoilNode
# phylo2node[phylo2node[,2] %in% SigPyOMTaxa$Name,]->SigPyOMNode
# SigSoilNode$Node->Node
# for(i in 1:dim(SigSoilNode)[1]){
#   p+geom_hilight(Node[i], fill=color3, alpha=0.5)->p
# }
# p

#an alternatinative for phylogenetic tree
# read.csv("phylodistance.csv",header = TRUE)->phylodist
# as.dist(phylodist)->phylodist
# cmdscale(phylodist,k=8,eig=TRUE)->phylo.MDS
# phylo.MDS$eig[1:20] #scree-plot shows k=5 is good
# library(ggplot2)
# phylo.MDS$points->phylopoint
# phylopoint<-as.data.frame(phylopoint)
# cbind(phylopoint,rownames(phylopoint))->phylopoint
# phylopoint<-phylopoint[order(phylopoint$`rownames(phylopoint)`),]
# cbind(phylopoint,anno.data)->phylopoint
# ggplot(data=phylopoint,aes(x=V1,y=V2))+geom_point(aes(color=response))
# rm(phylodist,phylo.MDS,phylopoint)#delete them just because not good

#testing code using PyNast
# ape::read.tree("16S.pynast.tre")->pynastTre
# ape::root(pynastTre, 1, r = TRUE)->root.pynastTre
# AnnotateTre(root.pynastTre,anno.data[anno.data$X.OTU.ID%in%root.pynastTre$tip.label,c(1,3)])->annotatedpynasttre
# attributes(annotatedpynasttre[[1]])$group->group
# as.character(group)->group
# unique(group)->taxa
# for(i in 1:length(taxa)) if(length(group[group==taxa[i]])<61) group[group==taxa[i]]<-"Other"
# as.factor(group)->attributes(annotatedpynasttre[[1]])$group
# ggtree(annotatedpynasttre[[1]],layout="fan",aes(color=group))+theme_minimal()+
#   theme(panel.grid = element_blank(),legend.position="right",
#         axis.text.x = element_blank(),axis.text.y = element_blank())->p.pynast
# #+xlim(0, 1.8)->p
# p.pynast %<+% anno.data + #geom_tippoint(aes(color=corr,subset=test), alpha=1,size=2) +
#   scale_color_manual(values = c(tol18rainbow[1:11],
#                                 "#D4D4D4",tol18rainbow[12:16],
#                                 "#888888",tol18rainbow[17:18]))




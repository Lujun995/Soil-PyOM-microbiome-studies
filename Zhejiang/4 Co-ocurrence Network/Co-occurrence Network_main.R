#data depends on D:\World Forest Data\Fuyang\data\.Rdata, ordination result and post-hoc selection result
library(ggplot2)
load("..\\.Rdata")
dget("NMDSresult")->w0.5NMDS.point
dget("anno.data")->anno.data
source("EasyNetwork V1.4.R")
source("MIC2p.R")
#partial regression with spearman correlation can control the type I error, even for an unbalanced design
#but it can only control data fitting model as X=f(Y)+Ci+e, where f(Y)is a linear or unlinear function of Y 
#and Ci is the main effect of different categories.
#asin-sqrt-transformation will partly solve the X=f(Y+Ci+e1)+e2 problem(variance depends on mean). 
#But there is no clear method yet.
asin(sqrt(rare/matrix(rep(rowSums(rare),ncol(rare)),nrow=nrow(rare),byrow=F)))->rare.transformed
filter.sparsity(rare.transformed)->rare.transformed.filtered
resid(lm(rare.transformed.filtered~Sample.Type))->rare.transformed.filtered.residual
cor(rare.transformed.filtered.residual,method="spearman")->rho
cor(rare.transformed.filtered.residual)->r
minerva::mine(rare.transformed.filtered.residual,n.cores=4)$MIC->MIC
r2p(r,n=43)->r.p
r2p(rho,n=43)->rho.p
MIC2p(MIC)->MIC.p
mat2vec(rho.p,triangle=T)->rho.p  #an approximate caculation
mat2vec(rho,triangle=T)->rho
mat2vec(MIC.p,triangle=T)->MIC.p
mat2vec(MIC,triangle=T)->MIC
mat2vec(r.p,triangle=T)->r.p
mat2vec(r,triangle=T)->r
p.adjust(r.p$value,method="BH")->adj.p
cbind(r.p,adj.p)->r.p
p.adjust(rho.p$value,method="BH")->adj.p
cbind(rho.p,adj.p)->rho.p
p.adjust(MIC.p$value,method="BH")->adj.p
cbind(MIC.p,adj.p)->MIC.p
#Edges and Nodes in the network
cbind(MIC.p,rho.p,rho)->edges
dim(edges[edges[,8]<0.001&edges[,4]<0.05,])[1]
edges[edges[,8]<0.001&edges[,4]<0.05,c(2,3,1,4,5,8,9)]->edges #FDR of spearman<0.001, that of MIC<0.05
names(edges)<-c("From","To","MIC.p","MIC.q","Spearman.p","Spearman.q","rho")
sign(edges$rho)->edges$rho
write.csv(edges,"edges.csv")
nodes<-c(as.matrix(edges[,c("From","To")]))
nodes<-nodes[!duplicated(nodes)]
anno.data[anno.data$X.OTU.ID %in% nodes,]->nodes
write.csv(nodes,"nodes.csv")
#relative abundance in different group
ggplot(w0.5NMDS.point)+geom_dotplot(aes(MDS1))+
  scale_color_manual(values = Sample.Type)+scale_fill_manual(values = Sample.Type)
group1<-which(w0.5NMDS.point$MDS1>0.2)
group2<-which(w0.5NMDS.point$MDS1>-0.1 & w0.5NMDS.point$MDS1<0.2)
group3<-which(w0.5NMDS.point$MDS1< -0.1)
variable<-as.character(nodes$X.OTU.ID)
tempfunc<-function(comm,group,variable){
  comm[group,variable]->temp
  apply(temp,MARGIN=2,mean)->MEAN
  apply((temp==0),MARGIN=2,sum)/length(group)->sparsity
  data.frame(mean=MEAN,sparsity=sparsity)
}
cbind(tempfunc(rare,group1,variable),tempfunc(rare,group2,variable),tempfunc(rare,group3,variable))->temp
write.csv(temp,"relativeabundance.csv")
#some visulazation
ggplot(w0.5NMDS.point)+
  geom_dotplot(aes(MDS1,color=Sample.Type,fill=Sample.Type),dotsize=0.5)+
  theme_classic(12) + theme(panel.grid = element_blank(),legend.position = "right",
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.line.y = element_blank()) + 
  labs(x="NMDS1")+
  scale_colour_manual(values = c(color1,color2,color3)) + 
  scale_fill_manual(values = c(color1,color2,color3))
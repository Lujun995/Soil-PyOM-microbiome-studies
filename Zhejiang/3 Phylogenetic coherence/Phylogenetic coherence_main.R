#import required functions, dataset and packages
library(ggtree)
library(ape)
library(picante)
library(caper)
library(ggplot2)
source('comparative.data.R')
source('phylo.d.R')
source('summary.cladehub.R')
source("SigClade4.R")
source("as.phylo.data.frame.R")
source('pSpecialists.R')
source("gheatmap.R")
source("AnnotateTre.R")
source("groupClade.R")
source("CountTaxa.R")
source("Eliminate.R")
source("getDescendants.R")
source("AnnotationNode.R")




#Coherent lineages (Table S2a)
load("..\\.Rdata")
dget("anno.data")->anno.data #requires data from specialist selection
anno.data$response->response
response[!anno.data$ifsignificant]<-0
response<-as.factor(response)
levels(response)<-c("Positive","NA","Negative")
anno.data$response<-response
anno.data<-anno.data[,-ncol(anno.data)]
response<-as.data.frame(response)
rownames(response)<-rownames(anno.data)
SigClade4(anno.data[,1:7],which(anno.data$response=="Positive"),which(anno.data$response=="Negative"),Bonferroni=0.05)->all.coherent
temp<-matrix(0, nrow=nrow(all.coherent),ncol=3)
colnames(temp)<-c("PyOM","Neutral","Soil")
cbind(all.coherent,temp)->all.coherent
for(i in 1:nrow(all.coherent)) {
  anno.data[,c(as.character(all.coherent$Phylogeny)[i],"response")]->temp
  summary(temp[temp[,1]==as.character(all.coherent$taxa)[i],2])[c("PyOM","NA's","Soil")]->all.coherent[i,5:7]
}
write.csv(all.coherent,"CoherenceUnits.csv")




#Visualization (Fig 7a)
load("..\\.Rdata")
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




#Fritz's 1-D (Fig 7a)
load("..\\.Rdata")
dget("anno.data")->anno.data
anno.data$response->response
response[!anno.data$ifsignificant]<-0
rep(0,length(response))->temp1->temp2
temp1[response==-1]<-1
temp2[response==1]<- 1
category<-data.frame(OTUID=anno.data$X.OTU.ID,
                     positives=temp1,negatives=temp2,response=response)
rm(temp1,temp2)

#phylo.D
summary.cladehub(annotatedtre[[2]])->LargeClade
D.Transformed<-rep(NA,nrow(LargeClade))
P<-rep(NA,nrow(LargeClade))
for(i in 1:length(LargeClade$Node))  try(
  {phylo.d(category,ape::extract.clade(node = LargeClade$Node[i],phy = annotatedtre[[1]]),
           "OTUID","positives",permut = 300)->D
    -D$DEstimate+1->D.Transformed[i]
    D$Pval1->P[i]
  }
)
positives=D.Transformed
cbind(LargeClade,positives,P)->LargeClade
rm(positives)
D.Transformed<-rep(NA,nrow(LargeClade))
P<-rep(NA,nrow(LargeClade))
for(i in 1:length(LargeClade$Node))  try(
  {phylo.d(category,ape::extract.clade(node = LargeClade$Node[i],phy = annotatedtre[[1]]),
           "OTUID","negatives",permut = 300)->D
    -D$DEstimate+1->D.Transformed[i]
    D$Pval1->P[i]
  }
)
negatives=D.Transformed
cbind(LargeClade,negatives,P)->LargeClade  #Actually, D in LargeClade is 1-D.original. 
#Negative->soil, positive->PyOM
rm(negatives)

phylo.d(category,phylotre,
        "OTUID","negatives",permut = 300)->Drecord  #1-D=0.5058
phylo.d(category,phylotre,
        "OTUID","positives",permut = 300)->Drecord2  #1-D=0.2633

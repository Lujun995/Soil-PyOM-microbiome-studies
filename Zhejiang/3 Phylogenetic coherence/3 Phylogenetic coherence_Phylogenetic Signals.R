#phylogenetic signal
library(picante)
library(caper)
source('comparative.data.R')
source('phylo.d.R')
source('summary.cladehub.R')

#also depends on result from SnowFlake.R
load("D:\\World Forest Data\\Fuyang\\data\\.Rdata")
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

#Bloomberg's K
# anno.data$response<-response
# anno.data<-anno.data[,-ncol(anno.data)]
# response<-as.data.frame(response)
# rownames(response)<-rownames(anno.data)
# 
# test2<-data.frame(response=as.numeric(response[test4$tip.label,]))
# rownames(test2)<-test4$tip.label
# Kcalc(test2$response,test4)



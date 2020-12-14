library(matrixStats)
source("DAA.R")
source("rho2p.R")
#It requires .RData in the parent folder
load("..\\.Rdata")

#Post-hoc selection with DAA
data.frame(meta=meta[c(-1,-3),])->meta.drop;rare[c(-1,-3),]->rare.drop
rare.drop[,-which(colSums(rare.drop)==0)]->rare.drop
permute_differential_analysis(meta.drop,t(rare.drop),'meta',perm.no = 9999)->result.drop
FDR.drop<-data.frame(raw=result.drop$p.raw,fdr=result.drop$p.adj.fdr,fwer=result.drop$p.adj.fwer)
rownames(FDR.drop)<-names(result.drop$p.raw)
FDR.drop[!is.na(FDR.drop$raw),]->FDR.drop
dim(FDR.drop[FDR.drop$fdr<0.01,]->b)

#Selected OTUs (not necessarily compatible with other files)
phylo[phylo$X.OTU.ID %in% rownames(FDR.drop),]->phylo.drop
rare.drop[,colnames(rare) %in% rownames(FDR.drop)]->rare.drop
if(sum(colnames(rare.drop)!=rownames(FDR.drop))+sum(phylo.drop$X.OTU.ID!=rownames(FDR.drop)))
  warning("Rare.drop, FDR.drop, phylo.drop don't match.")
which(phylo.drop$X.OTU.ID %in% rownames(b))->Sig.drop
which(!(phylo.drop$X.OTU.ID %in% rownames(b)))->Nul.drop
sign(apply(rare.drop,MARGIN = 2,cor,meta.drop,method="spearman"))->response
ifsignificant<-(phylo.drop$X.OTU.ID %in% rownames(b))
if(sum(names(response)!=phylo.drop$X.OTU.ID)==0)
  cbind(phylo.drop,response,ifsignificant)->anno.drop else
    warning("phylo.drop, response, ifsignificant don't match.")
write.csv(anno.drop[anno.drop$ifsignificant,],"SelectedOTUs.csv")

#Selected OTUs (compatible file (OTUs=6100), for downstream analysis) 
sign(apply(rare,MARGIN = 2,cor,meta,method="spearman"))->response
data.response<-rep(NA,dim(data)[2])
names(data.response)<-colnames(data)
data.response[names(response)]<-response
ifsignificant<-(phylo$X.OTU.ID %in% rownames(b))
if(sum(names(data.response)!=phylo$X.OTU.ID)==0)
  cbind(phylo,data.response,ifsignificant)->anno.data else
    warning("phylo, data.response, ifsignificant don't match.")
names(anno.data)[9]<-"response"
dput(anno.data,"..\\3 Phylogenetic coherence\\anno.data")
dput(anno.data,"..\\4 Co-ocurrence Network\\anno.data")

#description of selected OTUs
View(anno.data)
anno.data$response[!anno.data$ifsignificant]<-0
anno.data$response<-as.factor(anno.data$response)
levels(anno.data$response)<-c("Positive","Neutral","Negative") #Positive to PyOM
length(which(anno.data$ifsignificant)) # 1020 OTUs
length(which(anno.data$ifsignificant))/(dim(anno.data)[1]) #16.72% of total OTUs
sum(rare[,as.character(anno.data$X.OTU.ID[anno.data$ifsignificant])])/sum(rare) #56.62% of total abundance
length(which(anno.data$response=="Positive")) #557 OTUs
length(which(anno.data$response=="Negative")) #463 OTUs
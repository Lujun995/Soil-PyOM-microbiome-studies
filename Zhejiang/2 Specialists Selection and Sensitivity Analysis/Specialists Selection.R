library(matrixStats)
source("DAA.R")
source("rho2p.R")
load("D:\\World Forest Data\\Fuyang\\data\\.Rdata")
#Post-Hoc using DAA
# permute_differential_analysis(meta,t(data),'meta',perm.no = 9999)->result
# FDR<-data.frame(raw=result$p.raw,fdr=result$p.adj.fdr,fwer=result$p.adj.fwer)
# rownames(FDR)<-names(result$p.raw)
# dim(FDR[FDR$fdr<0.001,]->a)
# if(sum(colnames(data)!=rownames(FDR))+sum(phylo$X.OTU.ID!=rownames(FDR)))
#   warning("data, FDR, phylo don't match.")
# which(phylo$X.OTU.ID %in% rownames(a))->Sig
# which(!(phylo$X.OTU.ID %in% rownames(a)))->Nul
# ifsignificant<-(phylo$X.OTU.ID %in% rownames(a))
# sign(apply(data,MARGIN = 2,cor,meta,method="spearman"))->response
# if(sum(names(response)!=phylo$X.OTU.ID)==0)
#   cbind(phylo,response,ifsignificant)->anno  else
#     warning("response, phylo don't match.")
# write.csv(anno[anno$ifsignificant,],"original.csv")

#Post-hoc selection with DAA
data.frame(meta=meta[c(-1,-3),])->meta.drop;rare[c(-1,-3),]->rare.drop
rare.drop[,-which(colSums(rare.drop)==0)]->rare.drop
permute_differential_analysis(meta.drop,t(rare.drop),'meta',perm.no = 9999)->result.drop
FDR.drop<-data.frame(raw=result.drop$p.raw,fdr=result.drop$p.adj.fdr,fwer=result.drop$p.adj.fwer)
rownames(FDR.drop)<-names(result.drop$p.raw)
FDR.drop[!is.na(FDR.drop$raw),]->FDR.drop
dim(FDR.drop[FDR.drop$fdr<0.01,]->b)

#Selected OTUs (not necceray compatible with other files)
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

#compatible format (OTUs=6100) for downstream analysis 
sign(apply(rare,MARGIN = 2,cor,meta,method="spearman"))->response
data.response<-rep(NA,dim(data)[2])
names(data.response)<-colnames(data)
data.response[names(response)]<-response
ifsignificant<-(phylo$X.OTU.ID %in% rownames(b))
if(sum(names(data.response)!=phylo$X.OTU.ID)==0)
  cbind(phylo,data.response,ifsignificant)->anno.data else
    warning("phylo, data.response, ifsignificant don't match.")
names(anno.data)[9]<-"response"
dput(anno.data,"D:\\World Forest Data\\Fuyang\\data\\4PhylogenyAnalysis\\anno.data")
dput(anno.data,"D:\\World Forest Data\\Fuyang\\data\\6Network Analysis\\anno.data")
#the result is identical, hence no debug should be done
# sign(apply(rare,MARGIN = 2,cor,meta,method="spearman"))->temp1
# sign(apply(data,MARGIN = 2,cor,meta,method="spearman"))->temp2
# (phylo.rare$X.OTU.ID %in% rownames(b))->temp3
# (phylo$X.OTU.ID %in% rownames(b))->temp4
# temp1[temp3]->temp5
# temp2[temp4]->temp6
# sum((temp5-temp6)^2)

#description of selected OTUs
View(anno.data)
anno.data$response[!anno.data$ifsignificant]<-0
anno.data$response<-as.factor(anno.data$response)
levels(anno.data$response)<-c("Positive","Neutral","Negative") #Positive to PyOM
length(which(anno.data$ifsignificant)) # 1020 OTUs
length(which(anno.data$ifsignificant))/(dim(anno.data)[1]) #16.72% of total OTUs
#length(which(anno.drop$ifsignificant))/(dim(anno.drop)[1]) #17.85%
sum(rare[,as.character(anno.data$X.OTU.ID[anno.data$ifsignificant])])/sum(rare) #56.62% of total abundance
length(which(anno.data$response=="Positive")) #557 OTUs
#PyOM samples were more than soil samples, hence this comparison is meaningless
#sum(rare[,as.character(anno.data$X.OTU.ID[anno.data$response=="Positive"])])/sum(rare) #35.05% of total abundance
length(which(anno.data$response=="Negative")) #463 OTUs
#sum(rare[,as.character(anno.data$X.OTU.ID[anno.data$response=="Negative"])])/sum(rare) #21.56% of total abundance

#Other attemps
# weight=rep(1, times=length(meta$meta))
# weight[c(1,3)]<-c(0.01,0.01)
# permute_differential_analysis(meta,t(rare),'meta',perm.no = 9999, weights = weight)->result.weight
# FDR.weight<-data.frame(raw=result.weight$p.raw,fdr=result.weight$p.adj.fdr,fwer=result.weight$p.adj.fwer)
# dim(FDR.weight[FDR.weight$fdr<0.001,]->c)
# 
# rare.rank=apply(rare,MARGIN = 2, rank, ties.method="random")
# permute_differential_analysis(meta,t(rare.rank),'meta',perm.no = 9999,
#                               transform="none",size.factor =rep(1,45))->result.rank
# FDR.rank<-data.frame(raw=result.rank$p.raw,fdr=result.rank$p.adj.fdr,
#                      fwer=result.rank$p.adj.fwer)
# dim(FDR.rank[FDR.rank$fdr<0.001,]->d)
# 
# 
# apply(rare, MARGIN = 2, cor, method="spearman",y=meta$meta)->rho
# rho2p(rho,45)[,1]->p.raw
# BH.Holm.rank<-data.frame(p.raw=p.raw,fdr=p.adjust(p.raw,method="fdr"),fwer=p.adjust(p.raw,method="holm"))
# rm(p.raw,rho)
# dim(BH.Holm.rank[BH.Holm.rank$fdr<0.001,]->e)
# 
# venn.plot <- draw.quintuple.venn(
#   area1 = length(rownames(a)),
#   area2 = length(rownames(b)),
#   area3 = length(rownames(c)),
#   area4 = length(rownames(d)),
#   area5 = length(rownames(e)),
#   n12 = length(intersect(rownames(a),rownames(b))),
#   n13 = length(intersect(rownames(a),rownames(c))),
#   n14 = length(intersect(rownames(a),rownames(d))),
#   n15 = length(intersect(rownames(a),rownames(e))),
#   n23 = length(intersect(rownames(b),rownames(c))),
#   n24 = length(intersect(rownames(b),rownames(d))),
#   n25 = length(intersect(rownames(b),rownames(e))),
#   n34 = length(intersect(rownames(c),rownames(d))),
#   n35 = length(intersect(rownames(c),rownames(e))),
#   n45 = length(intersect(rownames(d),rownames(e))),
#   n123 = length(intersect(intersect(rownames(a),rownames(b)),rownames(c))),
#   n124 = length(intersect(intersect(rownames(a),rownames(b)),rownames(d))),
#   n125 = length(intersect(intersect(rownames(a),rownames(b)),rownames(e))),
#   n134 = length(intersect(intersect(rownames(a),rownames(c)),rownames(d))),
#   n135 = length(intersect(intersect(rownames(a),rownames(c)),rownames(e))),
#   n145 = length(intersect(intersect(rownames(a),rownames(d)),rownames(e))),
#   n234 = length(intersect(intersect(rownames(b),rownames(c)),rownames(d))),
#   n235 = length(intersect(intersect(rownames(b),rownames(c)),rownames(e))),
#   n245 = length(intersect(intersect(rownames(b),rownames(d)),rownames(e))),
#   n345 = length(intersect(intersect(rownames(c),rownames(d)),rownames(e))),
#   n1234 = length(intersect(intersect(intersect(rownames(a),rownames(b)),rownames(c)),rownames(d))),
#   n1235 = length(intersect(intersect(intersect(rownames(a),rownames(b)),rownames(c)),rownames(e))),
#   n1245 = length(intersect(intersect(intersect(rownames(a),rownames(b)),rownames(d)),rownames(e))),
#   n1345 = length(intersect(intersect(intersect(rownames(a),rownames(c)),rownames(d)),rownames(e))),
#   n2345 = length(intersect(intersect(intersect(rownames(b),rownames(c)),rownames(d)),rownames(e))),
#   n12345 = length(intersect(intersect(intersect(intersect(rownames(b),rownames(c)),rownames(d)),rownames(e)),rownames(a))),
#   category = c("ori", "drop", "weight", "rank", "BH"),
#   fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
#   cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
#   cat.cex = 2,
#   margin = 0.05,
#   cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 
#           1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
#   ind = TRUE)

#By additive
dim(FDR.drop[FDR.drop$fdr<1e-5,]) #253 OTUs
which(phylo.drop$X.OTU.ID %in% rownames(FDR.drop[FDR.drop$fdr<1e-5,]))->Sig.drop1
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=Sig.drop1,time=100)
dim(FDR.drop[FDR.drop$fdr<3e-5,]) #320 OTUs
which(phylo.drop$X.OTU.ID %in% rownames(FDR.drop[FDR.drop$fdr<3e-5,]))->Sig.drop2
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=Sig.drop2,time=100)
dim(FDR.drop[FDR.drop$fdr<1e-4,]) #413 OTUs
which(phylo.drop$X.OTU.ID %in% rownames(FDR.drop[FDR.drop$fdr<1e-4,]))->Sig.drop3
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=Sig.drop3,time=100)
dim(FDR.drop[FDR.drop$fdr<3e-4,]) #523 OTUs
which(phylo.drop$X.OTU.ID %in% rownames(FDR.drop[FDR.drop$fdr<3e-4,]))->Sig.drop4
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=Sig.drop4,time=100)
dim(FDR.drop[FDR.drop$fdr<1e-3,]) #628 OTUs
which(phylo.drop$X.OTU.ID %in% rownames(FDR.drop[FDR.drop$fdr<1e-3,]))->Sig.drop5
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=Sig.drop5,time=100)
dim(FDR.drop[FDR.drop$fdr<3e-3,]) #761 OTUs
which(phylo.drop$X.OTU.ID %in% rownames(FDR.drop[FDR.drop$fdr<3e-3,]))->Sig.drop6
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=Sig.drop6,time=100)
dim(FDR.drop[FDR.drop$fdr<1e-2,]) #1020 OTUs
which(phylo.drop$X.OTU.ID %in% rownames(FDR.drop[FDR.drop$fdr<1e-2,]))->Sig.drop7
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=Sig.drop7,time=100)
dim(FDR.drop[FDR.drop$fdr<3e-2,]) #1381
which(phylo.drop$X.OTU.ID %in% rownames(FDR.drop[FDR.drop$fdr<3e-2,]))->Sig.drop8
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=Sig.drop8,time=100)
dim(FDR.drop[FDR.drop$fdr<1e-1,]) #1954
which(phylo.drop$X.OTU.ID %in% rownames(FDR.drop[FDR.drop$fdr<1e-1,]))->Sig.drop9
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=Sig.drop9,time=100)
dim(FDR.drop[FDR.drop$fdr<3e-1,]) #3028
which(phylo.drop$X.OTU.ID %in% rownames(FDR.drop[FDR.drop$fdr<3e-1,]))->Sig.drop10
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=Sig.drop10,time=100)
dim(FDR.drop) #5714

#By differnence
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=setdiff(Sig.drop2,Sig.drop1),time=100)
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=setdiff(Sig.drop3,Sig.drop2),time=100)
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=setdiff(Sig.drop4,Sig.drop3),time=100)
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=setdiff(Sig.drop5,Sig.drop4),time=100)
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=setdiff(Sig.drop6,Sig.drop5),time=100)
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=setdiff(Sig.drop7,Sig.drop6),time=100)
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=setdiff(Sig.drop8,Sig.drop7),time=100)
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=setdiff(Sig.drop9,Sig.drop8),time=100)
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=setdiff(Sig.drop10,Sig.drop9),time=100)
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=setdiff(rownames(rare.drop),Sig.drop10),time=100)

library(picante)
library(ggplot2)
library(reshape2)
library(parallel)
#source("https://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
#library(devtools)
#install_github("Russel88/MicEco")
#library(MicEco)  #there is a bug in MicEco, don't use it
mat2vec<-function(x,triangle=FALSE,data.in=NULL){
  #to convert a matrix or data.frame to three vectors: value,rname,cname
  #especially when the matrix or data.frame is symmetric or triangle.
  rownames(x)->rnam
  colnames(x)->cnam
  if(!is.matrix(x)){
    if(is.data.frame(x)) as.matrix(x)->x
    else stop("Input is not a dataframe or matrix.")
  }
  #to test if it's triangle or sysmetric
  if(isSymmetric(x)){
    if(!identical(rnam,cnam)) warning("Row names is not equal 
                                      to Col names.")
    triangle=TRUE
    data.in="lower"
  }
  if(triangle){
    if(!identical(rnam,cnam)) warning("Row names is not equal 
                                      to Col names.")
    if(data.in=="lower") x[upper.tri(x,diag=TRUE)]<-NA
    else if(data.in=="upper") x[lower.tri(x,diag=TRUE)]<-NA
    else stop("You said it's a triangle matrix, but where is
              the data?")
  }
  #to convert	
  c(x)->value
  rep(rnam,time=length(cnam))->rowname
  c(matrix(rep(cnam,time=length(rnam)),byrow=T,nrow=length(rnam),
           ncol=length(cnam)))->colname
  data.frame(value,rowname,colname)->y
  if(triangle){
    y[!is.na(y$value),]->y
    rownames(y)<-1:length(rownames(y))
  }
  y
  }

#NTI
#abundance=T
ses.mntd(rare,cophenetic(root.phylotre),null.model ="taxa.labels",
         abundance.weighted=TRUE)->NTI.df
NTI<-data.frame(NTI=-NTI.df$mntd.obs.z,SampaleType=Sample.Type)
factor(NTI$SampaleType,levels = c("UnburntSoil","BurntSoil","PyOM"))->NTI$SampaleType
ggplot(data=NTI)+
  geom_boxplot(aes(y=NTI,x=SampaleType,fill=SampaleType))+
  scale_fill_manual(values=c(color1,color2,color3))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())
#abundance=F
# ses.mntd(rare,cophenetic(root.phylotre),null.model ="taxa.labels",abundance.weighted=F)->NTI.df2
# NTI2<-data.frame(NTI=-NTI.df2$mntd.obs.z,SampaleType=Sample.Type)
# factor(NTI2$SampaleType,levels = c("UnburntSoil","BurntSoil","PyOM"))->NTI2$SampaleType
# ggplot(data=NTI2)+
#   geom_boxplot(aes(y=NTI,x=SampaleType,fill=SampaleType))+
#   scale_fill_manual(values=c(color1,color2,color3))+
#   theme_bw(base_size = 12)+
#   theme(panel.grid = element_blank())
# #abundance=T, asin-root-transformed
# rare.t<-asin(sqrt(rare)/rowSums(rare)[1])
# ses.mntd(rare.t,cophenetic(root.phylotre),null.model ="taxa.labels",abundance.weighted=T)->NTI.df3
# NTI3<-data.frame(NTI=-NTI.df3$mntd.obs.z,SampleType=Sample.Type)
# factor(NTI3$SampleType,levels = c("UnburntSoil","BurntSoil","PyOM"))->NTI3$SampleType
# ggplot(data=NTI3)+
#   geom_boxplot(aes(y=NTI,x=SampleType,fill=SampleType))+
#   scale_fill_manual(values=c(color1,color2,color3))+
#   theme_bw(base_size = 12)+
#   theme(panel.grid = element_blank())
# dput(NTI3,"D:\\World Forest Data\\Fuyang\\data\\1Description\\alphaNTI.d")

#not betaNTI, really. Exactly, NRI
# betaNTI.df<-ses.mpd(rare,cophenetic(root.phylotre),
#                     null.model ="taxa.labels",abundance.weighted=TRUE)
# betaNTI<- -betaNTI.df$mpd.obs.z

#betaNTI
drop.tip(root.phylotre,root.phylotre$tip.label[!(root.phylotre$tip.label %in% colnames(rare))])->root.phylotre.d
cophenetic(root.phylotre.d)->root.phylotre.d.dist
ses.beta.mntd(rare,root.phylotre.d.dist,abundance.weighted = T)->betaNTI
# comdistnt(rare,root.phylotre.d.dist,abundance.weighted = T)->obs
# #not run
# cl <- makeCluster(detectCores()-2)  
# clusterEvalQ(cl,library(picante))
# clusterExport(cl,c("root.phylotre.d.dist","rare"))
# parSapply(cl, 1:100, function(i,...) { rare->temp; 
#   #  sample(colnames(rare.t))->colnames(temp);  #wrong, because I need to calculate pair-wisely, I should permute them individually.
#   t(apply(temp,MARGIN=1,sample))->temp
#   colnames(rare)->colnames(temp)
#   comdistnt(temp,root.phylotre.d.dist,abundance.weighted = T) } )->ran
# stopCluster(cl)
# apply(X = ran,MARGIN=1,mean)->ran.mean
# apply(X = ran,MARGIN=1,sd)->ran.sd
# betaNTI<-obs
# betaNTI[1:990]<- -(as.numeric(obs)-ran.mean)/ran.sd    #-betaMNTD, i.e. betaNTI
# as.matrix(betaNTI)->betaNTI
# diag(betaNTI)<-NA
# View(betaNTI)
median(abs(betaNTI[isPyOM,isPyOM]),na.rm = T)  #median=7.45
median(abs(betaNTI[isBurnt,isBurnt]),na.rm = T)  #median=10.83
median(abs(betaNTI[isUnburnt,isUnburnt]),na.rm = T)  #median=11.81
median(betaNTI[c(2,4:12),c(2,4:12)],na.rm = T) #remove the two outliers, 12.22
# rownames(betaNTI.asin.root)<-colnames(betaNTI.asin.root)<-as.character(Sample.Type)
# betaNTI.asin.root.m<-reshape2::melt(betaNTI.asin.root[c(isUnburnt,isBurnt,isPyOM),
#                         c(isUnburnt,isBurnt,isPyOM)])
# dput(betaNTI.asin.root.m,"D:\\World Forest Data\\Fuyang\\data\\1Description\\betaNTI.d")
# betaNTI.asin.root.m.r<-reshape2::melt(betaNTI.asin.root[c(2,4:12,isBurnt,isPyOM),
#                                                       c(2,4:12,isBurnt,isPyOM)])
# betaNTI.asin.root.m.plot<-
#   betaNTI.asin.root.m[betaNTI.asin.root.m$Var1==betaNTI.asin.root.m$Var2 & !is.na(betaNTI.asin.root.m$value),]
# betaNTI.asin.root.m.r.plot<-
#   betaNTI.asin.root.m.r[betaNTI.asin.root.m.r$Var1==betaNTI.asin.root.m.r$Var2 & !is.na(betaNTI.asin.root.m.r$value),]
# ggplot(data=betaNTI.asin.root.m.plot,aes(x=Var1,y=value))+geom_violin(aes(fill=Var1))+geom_boxplot(width=0.075)
# ggplot(data=betaNTI.asin.root.m.r.plot,aes(x=Var1,y=value))+geom_violin(aes(fill=Var1))+geom_boxplot(width=0.075)
# wilcox.test(betaNTI.asin.root.m.plot[betaNTI.asin.root.m.plot$Var1=="UnburntSoil",3],
#             betaNTI.asin.root.m.plot[betaNTI.asin.root.m.plot$Var1=="BurntSoil",3])
# wilcox.test(betaNTI.asin.root.m.r.plot[betaNTI.asin.root.m.r.plot$Var1=="UnburntSoil",3],
#             betaNTI.asin.root.m.r.plot[betaNTI.asin.root.m.r.plot$Var1=="BurntSoil",3])

mat2vec(betaNTI[isUnburnt,isUnburnt])->betaNTI.U
mat2vec(betaNTI[isPyOM,isPyOM])->betaNTI.P
mat2vec(betaNTI[isBurnt,isBurnt])->betaNTI.B
mat2vec(betaNTI[isUnburnt,isBurnt])->betaNTI.UB
mat2vec(betaNTI[isUnburnt,isPyOM])->betaNTI.UP
mat2vec(betaNTI[isBurnt,isPyOM])->betaNTI.BP
betaNTI.U[,2]<-"U"; betaNTI.P[,2]<-"P"; betaNTI.B[,2]<-"B"
betaNTI.UB[,2]<-"UB"; betaNTI.UP[,2]<-"UP"; betaNTI.BP[,2]<-"BP"
rbind(betaNTI.U[,1:2],betaNTI.B[,1:2],betaNTI.P[,1:2],
      betaNTI.UB[,1:2],betaNTI.UP[,1:2],betaNTI.BP[,1:2])->betaNTI.v
betaNTI.v$rowname<-as.factor(betaNTI.v$rowname)
betaNTI.v$rowname<-factor(betaNTI.v$rowname,levels=c("U","B","P","UB","UP","BP"))
dunn.test::dunn.test(betaNTI.v$value,betaNTI.v$rowname,method = "bh")
ggplot(data=betaNTI.v, aes(x=rowname,y=value))+
  #  geom_segment(x=1,y=100,xend=2,yend=100,size=0.5,inherit.aes = F)+geom_segment(x=1,y=150,xend=3,yend=150,size=0.5,inherit.aes = F)+
#  geom_violin(aes(fill=rowname),linetype=0)+
  geom_violin(aes(fill=rowname),linetype=0)+geom_boxplot(width=0.2,lwd=0.2,outlier.size = 1.5)+
  #  geom_blank(aes(y = y_min))+geom_blank(aes(y = y_max))+  #just to set the limits of facet
  theme_bw(base_size=12)+ guides(fill=FALSE)+
  theme(panel.grid = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  geom_segment(y=2,x=0.5,yend=2,xend=6.5,linetype="dashed")+
  scale_fill_brewer(palette = 6,type = "qual")+ylab(expression(beta*"NTI"))+ylim(0,17)


ggplot(betaNIT.m,aes(Var1,Var2))+ 
  geom_tile(aes(fill=value),color = "lightgray")+
  guides(fill=guide_colorbar("betaNTI"))+
  scale_fill_gradientn(colors=c("seagreen","khaki1","firebrick2"),guide="colorbar") +
  theme_minimal(15) + 
  theme(panel.grid = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.x = element_text(angle = 270, hjust = 0,vjust=0.3,color="black"),
        axis.text.y = element_text(color="black")
  )->heatmap.betaNTI
heatmap.betaNTI
pheatmap::pheatmap(-betaNTI)

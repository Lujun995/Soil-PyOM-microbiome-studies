#require .RData from 2 Specialists Selection and Sensitivity Analysis
library(GUniFrac)
library(vegan)
library(ggplot2)
library(grid)
library(gridExtra)
library(picante)
library(reshape2)
library(data.table)
source("myfunc.R")
load("..\\2 Specialists Selection and Sensitivity Analysis\\.RData")

Specialist<-rare[sample(isPyOM,1e2,replace =T),Sig.rare]+
  rare[sample(isUnburnt,1e2,replace =T),Sig.rare] #averaging both PyOM and Unburnt samples
Generalist<-rare[sample(c(isPyOM,isUnburnt),1e2,replace =T),Nul.rare]+
  rare[sample(c(isUnburnt,isPyOM),1e2,replace =T),Nul.rare] #generalist can be from either
mock2<-matrix(data=0,nrow=1e2,ncol=dim(rare)[2])
mock2[,Sig.rare]<-Specialist
mock2[,Nul.rare]<-Generalist
mock2<-GUniFrac::Rarefy(mock2,20000)[[1]]
rownames(mock2)<-paste("Mock",1:1e2,sep ="")
rbind(rare,mock2)->mock2

isMock=46:145
Sample.Type=factor(c(rep("UnburntSoil",12),rep("PyOM",21),rep("BurntSoil",12),rep("Mock",1e2)),
                   levels=c("UnburntSoil","BurntSoil","Mock","PyOM"))

vegan::diversity(mock2)->shannon2
vegan::specnumber(mock2)->OTUnum2

ses.mntd(mock2,cophenetic(root.phylotre),
         null.model ="taxa.labels",abundance.weighted=TRUE,runs = 999)->NTI.df2
NTI2<-data.frame(NTI=-NTI.df2$mntd.obs.z,SampaleType=Sample.Type)
NTI2$NTI[isBurnt]<- -NTI.df$mntd.obs.z[isBurnt] #NTI2$NTI[isBurnt] is identical to the one used before

ks.test(shannon2[isMock],shannon2[isBurnt])  #D=0.2400, p=0.50
ks.test(OTUnum2[isMock],OTUnum2[isBurnt])  #D=0.1767, p=0.89
ks.test(NTI2$NTI[isMock],NTI2$NTI[isBurnt])  #D=0.3133, p=0.20

GUniFrac(mock2,root.phylotre)$unifracs[,,"d_0.5"]->mock2.GU
as.data.frame(metaMDS(mock2.GU,k=3)$points)->mock2.NMDS #stress=0.14

ks.test(mock2.NMDS$MDS1[isMock],mock2.NMDS$MDS1[isBurnt])  #D=0.2833, p=0.30
ks.test(mock2.NMDS$MDS2[isMock],mock2.NMDS$MDS2[isBurnt])  #D=0.7067, p<0.0001
ks.test(mock2.NMDS$MDS3[isMock],mock2.NMDS$MDS3[isBurnt])  #D=0.2967, p=0.2504


myfunc12(mock2.NMDS,Sample.Type, c(color3,color1,"gray30",color2),-0.4,0.4,-0.4,0.4)->mock2pic1
myfunc13(mock2.NMDS,Sample.Type, c(color3,color1,"gray30",color2),-0.4,0.4,-0.3,0.3)->mock2pic2

comparison<-data.frame(SampleType=Sample.Type,
                       Shannon=shannon2,
                       OTUNum=OTUnum2,
                       aNTI=NTI2$NTI)
comparison[c(isBurnt,isMock),]->comparison
ggplot(data=comparison, aes(x=SampleType,y=Shannon))+
  geom_violin(aes(fill=SampleType),linetype=0)+geom_boxplot(width=0.2,lwd=0.2)+
  theme_bw()+ guides(fill=FALSE)+
  theme(panel.grid = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  ylab("Shannon-Wiener")+ylim(5.6,7.0)+
  scale_fill_manual(values=c(color1,"gray50"))->com1
ggplot(data=comparison, aes(x=SampleType,y=OTUNum))+
  geom_violin(aes(fill=SampleType),linetype=0)+geom_boxplot(width=0.2,lwd=0.2)+
  theme_bw()+ guides(fill=FALSE)+
  theme(panel.grid = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  ylab("Observed OTUs")+ylim(1400,3000)+
  scale_fill_manual(values=c(color1,"gray50"))->com2
ggplot(data=comparison, aes(x=SampleType,y=aNTI))+
  geom_violin(aes(fill=SampleType),linetype=0)+geom_boxplot(width=0.2,lwd=0.2)+
  theme_bw()+ guides(fill=FALSE)+
  theme(panel.grid = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  ylab(expression(alpha*"NTI"))+scale_y_continuous(labels=scaleFUN,limits=c(3.25,8.5))+
  scale_fill_manual(values=c(color1,"gray50"))->com3

plots<-list(com1,com2,com3,mock2pic1,mock2pic2)
g <- ggplotGrob(mockpic2 + theme(legend.position = "bottom"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$width)
layout<-matrix(c(1,2,3,4,4,5,5),nrow=1,byrow = T)
MockPic<-grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],
                      layout_matrix=layout,nrow=1)
combined<-arrangeGrob(MockPic,legend,ncol = 1,heights = unit.c(unit(1, "npc") - lheight, lheight))
grid.newpage()
#Fig 6 a&b
grid.draw(combined)

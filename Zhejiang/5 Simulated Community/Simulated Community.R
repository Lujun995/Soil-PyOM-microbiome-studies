library(GUniFrac)
library(vegan)
library(ggplot2)
library(grid)
library(gridExtra)
library(picante)
library(reshape2)
library(data.table)
source("myfunc.R")
#Determinisitic and stochastic dispersal 
mock1<-matrix(data=0,nrow=1e2,ncol=dim(rare)[2])
mock1<-rare[sample(isPyOM,1e2,replace =T),]+rare[sample(isUnburnt,1e2,replace =T),]
mock1<-GUniFrac::Rarefy(mock1,20000)[[1]]
rownames(mock1)<-paste("Mock",1:1e2,sep ="")
rbind(rare,mock1)->mock1

isMock=46:145
Sample.Type=factor(c(rep("UnburntSoil",12),rep("PyOM",21),rep("BurntSoil",12),rep("Mock",1e2)),
                   levels=c("UnburntSoil","BurntSoil","Mock","PyOM"))

vegan::diversity(mock1)->shannon1
vegan::specnumber(mock1)->OTUnum1
ses.mntd(mock1,cophenetic(root.phylotre),
         null.model ="taxa.labels",abundance.weighted=TRUE,runs = 999)->NTI.df1
data.frame(SampleType=Sample.Type,NTI=-NTI.df1$mntd.obs.z)->NTI1
NTI1$NTI[isBurnt]<- -NTI.df$mntd.obs.z[isBurnt] #noticeably, NTI1$NTI[isBurnt] is identical to the one used in Fig 3

ks.test(shannon1[isMock],shannon1[isBurnt])  #D=0.2567, p=0.41
ks.test(OTUnum1[isMock],OTUnum1[isBurnt])  #D=0.2933, p=0.32
ks.test(NTI1$NTI[isMock],NTI1$NTI[isBurnt])  #D=0.2933, p=0.26

GUniFrac(mock1,root.phylotre)$unifracs[,,"d_0.5"]->mock1.GU
as.data.frame(metaMDS(mock1.GU,k=3)$points)->mock1.NMDS #stress=0.12

ks.test(mock1.NMDS$MDS1[isMock],mock1.NMDS$MDS1[isBurnt])  #D=0.33, p=0.1548
t.test(mock1.NMDS$MDS1[isMock],mock1.NMDS$MDS1[isBurnt])  #p=0.4756
ks.test(mock1.NMDS$MDS2[isMock],mock1.NMDS$MDS2[isBurnt])  #D=0.57, p=0.0009
t.test(mock1.NMDS$MDS2[isMock],mock1.NMDS$MDS2[isBurnt])  #p=0.0008
ks.test(mock1.NMDS$MDS3[isMock],mock1.NMDS$MDS3[isBurnt])  #D=0.40, p=0.0444
t.test(mock1.NMDS$MDS3[isMock],mock1.NMDS$MDS3[isBurnt])  #p=0.0353

ggplot(mock1.NMDS, aes(y=MDS2, x=MDS1, colour = Sample.Type)) +
  geom_point(size = 2,aes(alpha=Sample.Type)) +
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_manual(name="SampleType",values=c(color3,color1,"gray50",color2)) +
  scale_alpha_manual(name="SampleType",values=c(0.3,1,1,0.3))+
  labs(y = "NMDS2", x = "NMDS1")->mock1pic1.legend
myfunc12(mock1.NMDS,Sample.Type, c(color3,color1,"gray50",color2),-0.4,0.4,-0.4,0.35)->mock1pic1
# ggplot(mock1.NMDS, aes(y=MDS3, x=MDS1, colour = Sample.Type)) +
#   geom_point(size = 2,aes(alpha=Sample.Type)) + 
#   theme_bw() + theme(panel.grid = element_blank()) + 
#   scale_alpha_manual(values=c(0.3,1,1,0.3))+
#   scale_colour_manual(values=c(color3,color1,"gray50",color2)) + 
#   labs(y = "NMDS3", x = "NMDS1")->mock1pic2
myfunc13(mock1.NMDS,Sample.Type, c(color3,color1,"gray50",color2),-0.4,0.4,-0.2,0.25)->mock1pic2

data.frame(sample=Sample.Type,NMDS1=mock1.NMDS$MDS1,Shannon=shannon1)->NMDS.vs.Shannon1
NMDS.vs.Shannon1[34:145,]->NMDS.vs.Shannon1
ggplot(data=NMDS.vs.Shannon1,aes(x=NMDS1,y=Shannon))+
  geom_point(aes(color=sample))+
  stat_smooth(method="lm", formula=y~x, size=0.75, aes(color=sample))+
  theme_bw()+
  scale_colour_manual(values=c(color1,"gray50"))+
  ylab("Shannon-Wiener")+
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        axis.text.y = element_text(color="black"))->mock1pic4


#Deterministic dispersal
#Post-hoc selection with DAA
data.frame(meta=meta[c(-1,-3),])->meta.drop;rare[c(-1,-3),]->rare.drop
rare.drop[,-which(colSums(rare.drop)==0)]->rare.drop
permute_differential_analysis(meta.drop,t(rare.drop),'meta',perm.no = 9999)->result.drop
FDR.drop<-data.frame(raw=result.drop$p.raw,fdr=result.drop$p.adj.fdr,fwer=result.drop$p.adj.fwer)
rownames(FDR.drop)<-names(result.drop$p.raw)
FDR.drop[!is.na(FDR.drop$raw),]->FDR.drop
dim(FDR.drop[FDR.drop$fdr<0.05,]->b)

Sig.rareII<-which(colnames(rare) %in% rownames(b))
Nul.rareII<-setdiff(1:(dim(rare)[2]),Sig.rareII)

Specialist<-rare[sample(isPyOM,1e2,replace =T),Sig.rareII]+rare[sample(isUnburnt,1e2,replace =T),Sig.rareII]
Generalist<-rare[sample(c(isPyOM,isUnburnt),1e2,replace =T),Nul.rareII]+
  rare[sample(c(isUnburnt,isPyOM),1e2,replace =T),Nul.rareII]
#Generalist<-rare[sample(c(isPyOM,isUnburnt),1e3,replace =T),Nul.rare]*2
mock2<-matrix(data=0,nrow=1e2,ncol=dim(rare)[2])
mock2[,Sig.rareII]<-Specialist
mock2[,Nul.rareII]<-Generalist
mock2<-GUniFrac::Rarefy(mock2,20000)[[1]]
#mock<-GUniFrac::Rarefy(rare[sample(isPyOM,120,replace =T),]+rare[sample(isUnburnt,120,replace =T),],20000)[[1]]
rownames(mock2)<-paste("Mock",1:1e2,sep ="")
rbind(rare,mock2)->mock2

isMock=46:145
Sample.Type=factor(c(rep("UnburntSoil",12),rep("PyOM",21),rep("BurntSoil",12),rep("Mock",1e2)),
                   levels=c("UnburntSoil","BurntSoil","Mock","PyOM"))

vegan::diversity(mock2)->shannon2
vegan::specnumber(mock2)->OTUnum2
#ses.mntd(rare,cophenetic(root.phylotre),
#         null.model ="taxa.labels",abundance.weighted=TRUE,runs = 999)->NTI.df 
ses.mntd(mock2,cophenetic(root.phylotre),
         null.model ="taxa.labels",abundance.weighted=TRUE,runs = 999)->NTI.df2
NTI2<-data.frame(NTI=-NTI.df2$mntd.obs.z,SampaleType=Sample.Type)
NTI2$NTI[isBurnt]<- -NTI.df$mntd.obs.z[isBurnt] #noticeably, NTI2$NTI[isBurnt] is identical to the one used in Fig 3

ks.test(shannon2[isMock],shannon2[isBurnt])  #D=0.2567, p=0.41
ks.test(OTUnum2[isMock],OTUnum2[isBurnt])  #D=0.1433, p=0.98
ks.test(NTI2$NTI[isMock],NTI2$NTI[isBurnt])  #D=0.4633, p=0.01

GUniFrac(mock2,root.phylotre)$unifracs[,,"d_0.5"]->mock2.GU
as.data.frame(metaMDS(mock2.GU,k=3)$points)->mock2.NMDS #stress=0.14

ks.test(mock2.NMDS$MDS1[isMock],mock2.NMDS$MDS1[isBurnt])  #D=0.23, p=0.5295
t.test(mock2.NMDS$MDS1[isMock],mock2.NMDS$MDS1[isBurnt])  #p=0.6359
ks.test(mock2.NMDS$MDS2[isMock],mock2.NMDS$MDS2[isBurnt])  #D=0.63, p=0.0001
t.test(mock2.NMDS$MDS2[isMock],mock2.NMDS$MDS2[isBurnt])  #p=0.0004
ks.test(mock2.NMDS$MDS3[isMock],mock2.NMDS$MDS3[isBurnt])  #D=0.46, p=0.0130
t.test(mock2.NMDS$MDS3[isMock],mock2.NMDS$MDS3[isBurnt])  #p=0.0129

# ggplot(mock2.NMDS, aes(y=MDS2, x=MDS1, colour = Sample.Type)) +
#   geom_point(size = 2,aes(alpha=Sample.Type)) + 
#   theme_bw() + theme(panel.grid = element_blank()) + 
#   scale_colour_manual(name="SampleType",values=c(color3,color1,"gray30",color2)) + 
#   scale_alpha_manual(name="SampleType",values=c(0.3,1,1,0.3))+
#   labs(y = "NMDS2", x = "NMDS1")->mock2pic1
myfunc12(mock2.NMDS,Sample.Type, c(color3,color1,"gray30",color2),-0.4,0.4,-0.3,0.4)->mock2pic1
# ggplot(mock2.NMDS, aes(y=MDS3, x=MDS1, colour = Sample.Type)) +
#   geom_point(size = 2,aes(alpha=Sample.Type)) + 
#   theme_bw() + theme(panel.grid = element_blank()) + 
#   scale_alpha_manual(values=c(0.3,1,1,0.3))+
#   scale_colour_manual(values=c(color3,color1,"gray30",color2)) + 
#   labs(y = "NMDS3", x = "NMDS1")->mock2pic2
myfunc13(mock2.NMDS,Sample.Type, c(color3,color1,"gray30",color2),-0.4,0.4,-0.2,0.3)->mock2pic2

data.frame(sample=Sample.Type,NMDS1=mock2.NMDS$MDS1,Shannon=shannon2)->NMDS.vs.Shannon2
NMDS.vs.Shannon2[34:145,]->NMDS.vs.Shannon2
ggplot(data=NMDS.vs.Shannon2,aes(x=NMDS1,y=Shannon))+
  geom_point(aes(color=sample))+
  stat_smooth(method="lm", formula=y~x, size=0.75, aes(color=sample))+
  theme_bw()+
  scale_colour_manual(values=c(color1,"gray30"))+
  ylab("Shannon-Wiener")+
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        axis.text.y = element_text(color="black"))->mock2pic4

{
  isMock1<-46:145;isMock2<-146:245
comparison<-data.frame(SampleType=c(Sample.Type,Sample.Type[isMock]),
                       shannon=c(shannon1,shannon2[isMock]),
                       OTUNum=c(OTUnum1,OTUnum2[isMock]),
                       aNTI=c(NTI1$NTI,NTI2$NTI[isMock]))
comparison$SampleType[isBurnt]<-rep("BurntSoil",length(isBurnt))
comparison$SampleType[isMock1]<-rep("Model DS",length(isMock1))
comparison$SampleType[isMock2]<-rep("Model D",length(isMock2))
comparison<-comparison[34:245,]

ggplot(data=comparison, aes(x=SampleType,y=shannon))+
  geom_violin(aes(fill=SampleType),linetype=0)+geom_boxplot(width=0.2,lwd=0.2)+
  theme_bw()+ guides(fill=FALSE)+
  theme(panel.grid = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  ylab("Shannon-Wiener")+ylim(5.6,7.0)+
  scale_fill_manual(values=c(color1,"gray30","gray50"))->com1
ggplot(data=comparison, aes(x=SampleType,y=OTUNum))+
  geom_violin(aes(fill=SampleType),linetype=0)+geom_boxplot(width=0.2,lwd=0.2)+
  theme_bw()+ guides(fill=FALSE)+
  theme(panel.grid = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  ylab("Observed OTUs")+ylim(1400,3000)+
  scale_fill_manual(values=c(color1,"gray30","gray50"))->com2
scaleFUN <- function(x) sprintf("%.1f", x)
ggplot(data=comparison, aes(x=SampleType,y=aNTI))+
  geom_violin(aes(fill=SampleType),linetype=0)+geom_boxplot(width=0.2,lwd=0.2)+
  theme_bw()+ guides(fill=FALSE)+
  theme(panel.grid = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  ylab(expression(alpha*"NTI"))+scale_y_continuous(labels=scaleFUN,limits=c(3.25,8.5))+
  scale_fill_manual(values=c(color1,"gray30","gray50"))->com3
}
{
position="bottom"
plots<-list(com1,com2,com3,
            mock2pic1,mock2pic2,mock2pic4,
            mock1pic1,mock1pic2,mock1pic4)
g <- ggplotGrob(mock1pic1.legend + theme(legend.position = position))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$width)
gl<-plots
#gl <- lapply(plots, function(x) x + theme(legend.position="none"))

MockPic<-grid.arrange(gl[[1]],gl[[2]],gl[[3]],gl[[4]],gl[[5]],gl[[6]],gl[[7]],gl[[8]],gl[[9]],
                      nrow=3)

combined <- switch(position,
                   "bottom" = arrangeGrob(MockPic,
                                          legend,
                                          ncol = 1,
                                          heights = unit.c(unit(1, "npc") - lheight, lheight)),
                   "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                         legend,
                                         ncol = 2,
                                         widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
grid.newpage()
grid.draw(combined)
}



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
myfunc12(mock2.NMDS,Sample.Type, c(color3,color1,"gray50",color2),-0.4,0.4,-0.3,0.4)->mock2pic1
myfunc13(mock2.NMDS,Sample.Type, c(color3,color1,"gray50",color2),-0.4,0.4,-0.2,0.3)->mock2pic2
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
grid.draw(combined)

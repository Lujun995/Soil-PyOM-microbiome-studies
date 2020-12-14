#dependent on .Rdata from the parent folder
load("..\\.Rdata")

#import required packages and functions
library(ggplot2)
library(reshape2)
library(dunn.test)
library(ape)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(vegan)
library(picante)
library(data.table)
library(GUniFrac)
library(vegan)
library(picante)
library(parallel)
source("grid_arrange_shared_legend.R")
source("stat_bag.R")
scaleFUN <- function(x) sprintf("%.1f", x)
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
myfunc2<-function(data,Sample.Type,color){
  laymat=matrix(c(rep(c(rep.int(1,3),NA),times=1),
                  rep(c(rep.int(2,3),4),times=3),
                  rep(c(rep.int(3,3),5),times=3),
                  rep(6,8)),
                nrow=4,byrow = FALSE)
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  ggplot(data,aes(MDS1,colour=Sample.Type,fill=Sample.Type))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw(12) + theme(panel.grid = element_blank(),
                         axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank(),
                         axis.text.y=element_text(angle = 90,hjust = 0.5),
                         legend.position="none")+
    scale_y_continuous(labels=scaleFUN)+labs(x = "NMDS1")+
    scale_fill_manual(values = color)+scale_color_manual(values = color)+coord_flip()+
    xlim(-0.4,0.4)->pic1
  
  ggplot(data, aes(MDS2, MDS1, colour = Sample.Type, fill = Sample.Type)) +
    geom_point(size = 2) + 
    stat_bag(prop = 1, alpha = 0.4, show.legend = T) + 
    theme_bw(12) + theme(panel.grid = element_blank(),
                         #                         legend.key.height=unit(0.35,"inch"),
                         axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank()) + 
    scale_color_manual(values = color) + 
    scale_fill_manual(values = color) + 
    ylim(-0.4,0.4)+xlim(-0.35,0.35)->pic2
  ggplot(data, aes(MDS3, MDS1, colour = Sample.Type, fill = Sample.Type)) +
    geom_point(size = 2) + 
    stat_bag(prop = 1, alpha = 0.4, show.legend = T) + 
    theme_bw(12) + theme(panel.grid = element_blank(), 
                         axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         legend.position="none") + 
    scale_color_manual(values = color) + 
    scale_fill_manual(values = color) +
    ylim(-0.4,0.4)+xlim(-0.25,0.2)->pic3
  
  ggplot(data,aes(MDS2,colour=Sample.Type,fill=Sample.Type))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw(12) + 
    theme(panel.grid = element_blank(), 
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="none")+
    scale_y_continuous(labels=scaleFUN)+
    scale_fill_manual(values = color)+scale_color_manual(values = color)+
    xlim(-0.35,0.35)+labs(x = "NMDS2")->pic4
  ggplot(data,aes(MDS3,colour=Sample.Type,fill=Sample.Type))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw(12) + theme(panel.grid = element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank(),
                         legend.position="none")+
    scale_y_continuous(labels=scaleFUN)+
    scale_fill_manual(values = color)+scale_color_manual(values = color)+
    xlim(-0.25,0.2)+labs(x = "NMDS3")->pic5
  
  g<-ggplotGrob(pic2)$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  pic2<-pic2+theme(legend.position="none")
  grid.arrange(pic1,pic2,pic3,pic4,pic5,legend,layout_matrix=laymat)
}
myfunc3<-function(data,Sample.Type,color){
  laymat=matrix(
        c(rep(c(4,2,2,2),3),
          rep(c(5,3,3,3),3),
          rep(c(NA,1,1,1),1)),
    nrow=7,byrow = TRUE)
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  ggplot(data,aes(MDS1,colour=Sample.Type,fill=Sample.Type))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw(12) + theme(panel.grid = element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank(),
                         legend.position="none")+
    scale_y_continuous(labels=scaleFUN)+
    scale_fill_manual(values = color)+scale_color_manual(values = color)+
    xlim(-0.4,0.4)+labs(x = "NMDS1")->pic1
  
  ggplot(data, aes(MDS1, MDS2, colour = Sample.Type, fill = Sample.Type)) +
    geom_point(size = 2) + 
    stat_bag(prop = 1, alpha = 0.4, show.legend = T) + 
    theme_bw(12) + theme(panel.grid = element_blank(),
                         #legend.key.height=unit(0.35,"inch"),
                         axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         legend.position="none") + 
    scale_color_manual(values = color) + 
    scale_fill_manual(values = color) + 
    xlim(-0.4,0.4)+ylim(-0.35,0.35)->pic2
  ggplot(data, aes(MDS1, MDS3, colour = Sample.Type, fill = Sample.Type)) +
    geom_point(size = 2) + 
    stat_bag(prop = 1, alpha = 0.4, show.legend = T) + 
    theme_bw(12) + theme(panel.grid = element_blank(), 
                         axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         legend.position=c(0.5,0.1),
                         legend.title=element_blank(),legend.direction= "horizontal") + 
    scale_color_manual(values = color) + 
    scale_fill_manual(values = color) +
    xlim(-0.4,0.4)+ylim(-0.3,0.2)->pic3
  
  ggplot(data,aes(MDS2,colour=Sample.Type,fill=Sample.Type))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw(12) + theme(panel.grid = element_blank(),
                         axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank(),
                         axis.text.y=element_text(angle = 90,hjust = 0.5),
                         legend.position="none")+
    scale_y_continuous(labels=scaleFUN)+labs(x = "NMDS2")+
    scale_fill_manual(values = color)+scale_color_manual(values = color)+coord_flip()+
    xlim(-0.35,0.35)->pic4
  ggplot(data,aes(MDS3,colour=Sample.Type,fill=Sample.Type))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw(12) + theme(panel.grid = element_blank(),
                         axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank(),
                         axis.text.y=element_text(angle = 90,hjust = 0.5),
                         legend.position="none")+
    scale_y_continuous(labels=scaleFUN)+labs(x = "NMDS3")+
    scale_fill_manual(values = color)+scale_color_manual(values = color)+coord_flip()+
    xlim(-0.3,0.2)->pic5
  
#  g<-ggplotGrob(pic2)$grobs
#  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#  pic2<-pic2+theme(legend.position="none")
  grid.arrange(pic1,pic2,pic3,pic4,pic5,layout_matrix=laymat)
}
isUnburnt<-1:12
isPyOM<-13:33
isBurnt<-34:45
factor(Sample.Type,levels = c("UnburntSoil","BurntSoil","PyOM"))->Sample.Type




#data input
diversity(rare)->shannon
shannon
kruskal.test(shannon~Sample.Type)
pd(rare,root.phylotre)$PD->faith
faith
pielou=shannon/log(vegan::specnumber(rare)) #Pielou's J for evenness
pielou
ObservedOTU=vegan::specnumber(rare)
D=vegan::specnumber(rare)/sqrt(rowSums(rare)) #Species number can be used as richness. Menhinick's index
ses.mntd(rare,cophenetic(root.phylotre),null.model ="taxa.labels",
         abundance.weighted=TRUE)->NTI.df #alphaNTI, abundance=T
NTI<-data.frame(NTI=-NTI.df$mntd.obs.z,SampaleType=Sample.Type)
factor(NTI$SampaleType,levels = c("UnburntSoil","BurntSoil","PyOM"))->NTI$SampaleType
NTI->aNTI
pH=c(5.112,4.422,5.188,4.160,4.187,4.025,4.028,4.081,3.895,4.009,3.844,4.011,
     NA,NA,NA,NA,NA,NA,NA,NA,NA,
     4.898,5.749,5.666,6.415,5.583,5.353,5.814,6.008,6.107,5.955,NA,6.186,
     4.729,4.670,4.720,5.088,5.555,4.694,5.367,4.676,4.680,4.678,4.782,4.711)
watercontent=c(17.34,15.42,17.18,34.73,38.98,40.98,35.11,38.70,42.71,33.37,32.49,37.53,
               NA,NA,NA,NA,NA,NA,NA,NA,NA,
               87.28,85.90,87.74,82.73,80.80,83.28,87.29,86.34,85.92,87.40,92.86,86.07,
               23.16,27.05,25.03,17.62,21.18,21.69,18.57,29.06,26.09,24.40,24.10,26.03)
DOC=c(4.85,11.4,1.60,39.28,50.45,65.97,84.02,41.37,68.56 ,50.52,93.45,74.20,
      rep(NA,21),
  17.85,24.65,17.3,9.95,10.15,14.75,15.65,33.15,21.65,8.25,24.6,21.95)#significant
TOC=c(13.40,55.76,10.50,85.13,104.95,143.30,89.89,101.75,171.22,112.90,111.42,108.90,
      rep(NA,21),
      78.25,65.99,109.10,92.76,55.78,67.60,67.69,76.76,38.19,39.13,115.40,141.70)#no difference




#Physichemical Properties and alpha diversity (Fig 2 a-f)
data.frame(sample=Sample.Type,pH=pH,Shannon.Wiener=shannon,Faith=faith,Pielou=pielou,
           ObservedOTU=ObservedOTU,alphaNTI=aNTI$NTI)->description1
melt(description1,id=1)->description1.melt
factor(description1$sample,levels = c("UnburntSoil","BurntSoil","PyOM"))->description1$sample
levels(description1.melt$sample)<-c("Unburnt\nSoil","Burnt\nSoil","PyOM")
data.table(description1.melt)->description1.melt
description1.melt[variable=="alphaNTI",y_max := 9]
description1.melt[variable=="Shannon.Wiener",y_max := 7]
description1.melt[variable=="pH",y_max := 7]
description1.melt[variable=="Faith",y_max := 200]
description1.melt[variable=="Pielou",y_max := 0.9]
description1.melt[variable=="ObservedOTU",y_max := 2500]
levels(description1.melt$variable)<-c("pH",expression("Shannon"*"-"*"Wiener"),expression("Faith's"~"Diversity"),
                                      expression("Pielou's"~"J"),"Observed~OTUs",expression(alpha*"NTI"))
ggplot(data=description1.melt, aes(x=sample,y=value))+
  geom_violin(aes(fill=sample),linetype=0)+geom_boxplot(width=0.2,lwd=0.2,outlier.size = 0.7)+
  geom_blank(aes(y = y_max))+  
  facet_wrap(~ variable, scales="free", nrow=1, labeller = label_parsed)+
  theme_bw(base_size=10)+ guides(fill=FALSE)+
  theme(panel.grid = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  scale_fill_manual(values=c(color3,color1,color2))->pic1

#SOC & DOC (Fig S6a)
melt(descriptionS1,id=1)->descriptionS1.melt
factor(descriptionS1.melt$sample,levels = c("UnburntSoil","BurntSoil","PyOM"))->descriptionS1.melt$sample
levels(descriptionS1.melt$sample)<-c("Unburnt\nSoil","Burnt\nSoil","PyOM")
data.table(descriptionS1.melt)->descriptionS1.melt
levels(descriptionS1.melt$variable)<-c(expression("Faith's"~"Diversity"),expression("Pielou's"~"J"),"Observed~OTUs",
                                     "TOC~(gC/kgDw)","DOC~(mgC/kgDw)")
ggplot(data=descriptionS1.melt, aes(x=sample,y=value))+
  #  geom_segment(x=1,y=100,xend=2,yend=100,size=0.5,inherit.aes = F)+geom_segment(x=1,y=150,xend=3,yend=150,size=0.5,inherit.aes = F)+
  geom_violin(aes(fill=sample),linetype=0)+geom_boxplot(width=0.2,lwd=0.2)+
  geom_blank(aes(y = y_min))+geom_blank(aes(y = y_max))+  #just to set the limits of facet
  facet_wrap(~ variable, scales="free", nrow=1, labeller = label_parsed)+
  theme_bw(base_size=12)+ guides(fill=FALSE)+
  theme(panel.grid = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  scale_fill_manual(values=c(color3,color1,color2))->picS1

#normality test failed hence wilcox-bonferroni was used
nortest::lillie.test(shannon[isUnburnt])
nortest::lillie.test(faith[isUnburnt])
nortest::lillie.test(pH[isUnburnt])
nortest::lillie.test(pielou[isBurnt])
nortest::lillie.test(ObservedOTU[isUnburnt])
dunn.test::dunn.test(description$pH,description$sample,"bh") #raw p-value should multiply by two("two-sided")
dunn.test::dunn.test(description$Shannon.Wiener,description$sample,"bh")
dunn.test::dunn.test(description$Faith,description$sample,"bh")
dunn.test::dunn.test(pielou,Sample.Type,"bh") 
dunn.test::dunn.test(ObservedOTU,Sample.Type,"bh")
pairwise.t.test(aNTI$NTI,Sample.Type,"fdr",pool.sd = FALSE)

#normality holds so t-test is OK
t.test(description$SoilOrganicCarbon~Sample.Type,na.action = "na.omit",var.equal = TRUE)
t.test(description$DissolvableOrganicCarbon~Sample.Type,na.action = "na.omit",var.equal = FALSE)




#PhylumAbundance (Fig 5a)
read.csv("PhylumAbundance.csv",header =TRUE)->phylumAbundance
cbind(Sample.Type,phylumAbundance)->phylumAbundance
phylumAbundance.m<-melt(phylumAbundance,id=c(1,2))
names(phylumAbundance.m)[c(2,3)]<-c("sample","Phylum")
factor(phylumAbundance.m$Sample.Type,levels = c("UnburntSoil","BurntSoil","PyOM"))->
  phylumAbundance.m$Sample.Type
phylumAbundance.m$value<-phylumAbundance.m$value/200
ggplot(phylumAbundance.m,aes(sample,value))+
  geom_col(aes(fill=Phylum),position = position_stack(reverse = TRUE))+
  facet_grid(~ Sample.Type, scales = "free_x", space = "free_x")+
  scale_fill_brewer(type="qual",palette = 3)+
  ylab("Relative Abundance (%)")+
  theme_bw(12)+ 
  theme(panel.grid = element_blank(), axis.title.x=element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"))->phylumabundance.pic




#Ordination (Fig 4a, Fig S3)
GUniFrac(rare,root.phylotre)->unifrac
unifrac$unifracs[,,"d_0.5"]->w0.5unifrac
#adonis(w0.5unifrac~Class1,permutations = 999999) #p<e-6,based on factor model
adonis(w0.5unifrac~meta$meta,permutations = 999999) #p<e-6,based on linear model, more precise
adonis(w0.5unifrac[isPyOM,isPyOM]~Class2) #p=0.894, Washed v.s. unwashed

metaMDS(comm = w0.5unifrac,k =3)->NMDS0.5uf #stress=0.06
as.data.frame(NMDS0.5uf$points)->NMDS0.5uf.point
NMDS0.5uf.point$MDS1<- -NMDS0.5uf.point$MDS1 #-MDS1 was used for convenience, the same direct as PyOM gradient
#Fig 4a
myfunc3(NMDS0.5uf.point,Sample.Type,color=c(color3,color1,color2))->ordinationAll
#Fig S3
myfunc2(NMDS0.5uf.point[isPyOM,],PyOM.Type,color=RColorBrewer::brewer.pal(10, "Set3")[c(6,10)])->ordinationSub 




#Correlations between NMDS1, alpha diversity and pH (Fig 3 left part)
Shannon.NMDS1<-data.frame(sample=Sample.Type,Shannon=shannon,NMDS1=NMDS0.5uf.point$MDS1)
summary(lm(shannon ~ poly(NMDS0.5uf.point$MDS1,degree=2)))
summary(lm(faith ~ poly(NMDS0.5uf.point$MDS1,degree=2)))
summary(lm(pielou ~ poly(NMDS0.5uf.point$MDS1,degree=2)))
summary(lm(ObservedOTU ~ poly(NMDS0.5uf.point$MDS1,degree=2)))
ggplot(data=Shannon.NMDS1,aes(x=NMDS1,y=Shannon))+
  geom_point(aes(color=sample))+
  stat_smooth(method="lm", formula=y~poly(x,2), size=0.75, color="gray15")+
  theme_bw(12)+
  scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  ylab("Shannon-Wiener")+
  theme(panel.grid = element_blank(), 
        legend.position = c(0.5,0.1),legend.direction= "horizontal",
        legend.background = element_rect(color = "black",size = 0.25, linetype = "solid"),
        axis.text.y = element_text(color="black"))->Shannon.vs.NMDS1

index<-which(!is.na(pH))
summary(lm(shannon[index] ~ poly(pH[index],degree=2)))  
qplot(pH,shannon,color=Sample.Type)+
  stat_smooth(method="lm", formula=y~poly(x,2), size=0.75, color="gray15")+
  theme_bw(12)+scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  theme(panel.grid = element_blank(), 
                     axis.text.y = element_text(color="black"))+
  labs(x="pH",y="Shannon-Wiener Index")->pHvsShannon

cor.test(pH,NMDS0.5uf.point$MDS1) #r0.82, p=1e-9
qplot(pH,NMDS0.5uf.point$MDS1,color=Sample.Type)+
  stat_smooth(method="lm", formula=y~poly(x,1), size=0.75, color="gray15")+
  theme_bw(12)+theme(panel.grid = element_blank(), 
                     axis.text.y = element_text(color="black"))+
  labs(x="pH",y="NMDS1")+scale_color_manual(name= "SampleType",values=c(color3,color1,color2))->NMDS1vspH

summary(lm(shannon ~ poly(NMDS0.5uf.point$MDS1,degree=2)))  
qplot(NMDS0.5uf.point$MDS1,shannon,color=Sample.Type)+
  stat_smooth(method="lm", formula=y~poly(x,2), size=0.75, color="gray15")+
  theme_bw(12)+scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  theme(panel.grid = element_blank(), 
        axis.text.y = element_text(color="black"))+
  labs(x="NMDS1",y="Shannon-Wiener Index")->MDS1vsShannon

summary(lm(pielou ~ poly(NMDS0.5uf.point$MDS1,degree=2)))
qplot(NMDS0.5uf.point$MDS1,pielou,color=Sample.Type)+
  stat_smooth(method="lm", formula=y~poly(x,2), size=0.75, color="gray15")+
  theme_bw(12)+scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  theme(panel.grid = element_blank(), 
        axis.text.y = element_text(color="black"))+
  labs(x="NMDS1",y="Pielou's J")->MDS1vsPielou

summary(lm(ObservedOTU ~ poly(NMDS0.5uf.point$MDS1,degree=2))) 
qplot(NMDS0.5uf.point$MDS1,ObservedOTU,color=Sample.Type)+
  stat_smooth(method="lm", formula=y~poly(x,2), size=0.75, color="gray15")+
  theme_bw(12)+scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  theme(panel.grid = element_blank(), 
        axis.text.y = element_text(color="black"))+
  labs(x="NMDS1",y="Observed OTUs")->MDS1vsOTUnum

summary(lm(faith ~ poly(NMDS0.5uf.point$MDS1,degree=2))) 
qplot(NMDS0.5uf.point$MDS1,faith,color=Sample.Type)+
  stat_smooth(method="lm", formula=y~poly(x,2), size=0.75, color="gray15")+
  theme_bw(12)+scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  theme(panel.grid = element_blank(), 
        axis.text.y = element_text(color="black"))+
  labs(x="NMDS1",y="Faith's Diversity")->MDS1vsFaith
#Fig 3 (left part)
grid_arrange_shared_legend(NMDS1vspH,pHvsShannon,MDS1vsShannon,MDS1vsFaith,MDS1vsPielou,
                           MDS1vsOTUnum,ncol=2,nrow = 3,position = "bottom")




#Procrustes test (Fig 4c)
vegan::protest(X=NMDS0.5uf.point[22:33,],Y=NMDS0.5uf.point[34:45,],
               permutation=99999,symmetric=FALSE)->protest1 #p=0.003, R=0.68

qplot(x=NMDS0.5uf.point[1:45,1],y=NMDS0.5uf.point[1:45,3],aes(color=Sample.Type),size=I(2))+
  geom_segment(aes(xend=NMDS0.5uf.point[22:33,1],yend=NMDS0.5uf.point[22:33,3],
                   x=NMDS0.5uf.point[34:45,1],y=NMDS0.5uf.point[34:45,3]),arrow=arrow(length=unit(0.2,"cm")))+
  geom_bag(aes(x=NMDS0.5uf.point[1:12,1],y=NMDS0.5uf.point[1:12,3]),prop=1,alpha=0.4,color=color3,fill=color3)+
  theme_bw(12) + theme(panel.grid = element_blank()) + xlab("NMDS1")+ylab("NMDS3")+
  scale_color_manual(values = c(color3,color1,color2))

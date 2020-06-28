library(ggplot2)
library(reshape2)
library(dunn.test)
library(ape)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(vegan)
library(picante)
library(caper)
library(data.table)
library(GUniFrac)
library(vegan)
library(picante)
library(parallel)
library(matrixStats)
library(ggtree)
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\1&2 Description and Ordinations combination\\grid_arrange_shared_legend.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\1&2 Description and Ordinations combination\\stat_bag.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\3Post-hoc Selection &Sensitive Analysis\\DAA.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\3Post-hoc Selection &Sensitive Analysis\\rho2p.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\3Post-hoc Selection &Sensitive Analysis\\OTU.permutate.R")
source('C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\4PhylogenyAnalysis\\comparative.data.R')
source('C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\4PhylogenyAnalysis\\phylo.d.R')
source('C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\4PhylogenyAnalysis\\summary.cladehub.R')
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\4PhylogenyAnalysis\\SigClade3.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\4PhylogenyAnalysis\\as.phylo.data.frame.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\4PhylogenyAnalysis\\gheatmap.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\4PhylogenyAnalysis\\AnnotateTre.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\4PhylogenyAnalysis\\groupClade.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\4PhylogenyAnalysis\\CountTaxa.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\4PhylogenyAnalysis\\Eliminate.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\4PhylogenyAnalysis\\getDescendants.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\4PhylogenyAnalysis\\AnnotationNode.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\6Network Analysis\\EasyNetwork V1.4.R")
source("C:\\Users\\zlj\\Desktop\\Manuscript data\\Fuyang\\data\\10Simulations\\myfunc.R")
source('C:/Users/zlj/Desktop/Manuscript data/Fuyang/data/4PhylogenyAnalysis/SigClade4.R')
source('C:/Users/zlj/Desktop/Manuscript data/Fuyang/data/4PhylogenyAnalysis/pSpecialists.R')


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
    xlim(-0.3,0.45)+labs(x = "NMDS1")->pic1
  
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
    xlim(-0.3,0.45)+ylim(-0.35,0.35)->pic2
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
    xlim(-0.3,0.45)+ylim(-0.3,0.2)->pic3
  
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
FISHER<-function(sig,total,totalsig,all){
  tested<-matrix(c(sig,total-sig,totalsig-sig,all-totalsig-total+sig),ncol = 2,byrow = TRUE)
  fisher.test(tested,alternative = "greater")$p.value
}
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", 
                "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44",
                "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

#data input
read.csv("result.csv", header=TRUE)->data
data[,1]->rownames(data)
data[,-1]->data
as.matrix(data)->data
t(data)->data
min(rowSums(data))
GUniFrac::Rarefy(data, depth = 40000)$otu.tab.rff->rare
rare[,colSums(rare)!=0]->rare

read.csv("rep_tax_assignments.csv", header=TRUE)->anno
anno[anno$OTUID %in% colnames(rare),]->anno.rare
anno.rare[order(anno.rare$OTUID),]->anno.rare

ape::read.tree("16S.tre")->phylotre
ape::root(phylotre, 1, r = TRUE)->root.phylotre
ape::drop.tip(phy = root.phylotre, tip = root.phylotre$tip.label[!(root.phylotre$tip.label%in%anno.rare$OTUID)])->root.phylotre.rare

SampleType<-as.factor(rep(c("BurntSoil","PyOM","UnburntSoil"),each=12))
factor(SampleType,levels = c("UnburntSoil","BurntSoil","PyOM"))->SampleType
isBurnt<-1:12
isPyOM<-13:24
isUnburnt<-25:36
color1="#EE6A50"
color2="#00CD66"
color3="#5CACEE"

pH=c(4.56,4.98,4.93,4.89,4.9,4.79,4.84,5,5.03,4.66,
     5.17,5.11,4.52,5.51,4.88,4.69,5.17,4.79,5.46,5.45,
     4.82,5.08,5.58,5.34,3.99,4.06,3.93,4.27,3.99,4.08,
     4.16,4.09,4.06,4.28,4.2,4.19)
DOC=c(34.02,34.08,16.96,26.15,30.6,39.59,39.95,38.2,30.22,
      34.64,29.39,41.49,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,56.71,49.06,
      38.08,53.93,64.97,101.76,41.08,40.02,51.88,47.77,37.02,35.58)*5
watercontent=c(0.2376,0.1996 ,0.1493 ,0.1896,0.1038,0.1569,0.1397,0.1532,
         0.2140,0.1639 ,0.2571,0.1196,0.4526,0.6629,0.351,0.4666,0.4821,
         0.5476,0.4451,0.5299,0.3485,0.3505,0.4106,0.3805,0.2566,0.219,
         0.2152,0.228,0.2265,0.234,0.2309,0.2521,0.2214,0.2774,0.2411,
         0.2797)
TOC= c(38.59,47.28,47.72,39.76,47.11,64.97,38.91,40.6,
       42.45,36.49,69.84,46.95,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,29.14,
       41.33,22.51,19.09,22.3,49.24,27.43,33.21,26.99,33.5,28.42,36.65)

#description
diversity(rare)->shannon
shannon
pd(rare,root.phylotre)$PD->faith
faith
pielou=shannon/log(vegan::specnumber(rare)) #Pielou's J for evenness
pielou
ObservedOTU=vegan::specnumber(rare)
ObservedOTU
kruskal.test(ObservedOTU~SampleType)
ses.mntd(rare,cophenetic(root.phylotre),null.model ="taxa.labels",
         abundance.weighted=F)->NTI.df #alphaNTI, abundance=T
aNTI<-data.frame(NTI=-NTI.df$mntd.obs.z,SampleType=SampleType)
kruskal.test(aNTI$NTI~SampleType)
# data.frame(sample=SampleType,Shannon.Wiener=shannon,alphaNTI=aNTI$NTI,Faith=faith,ObservedOTU=ObservedOTU)->description1
# melt(description1,id=1)->description1.melt
# factor(description1$sample,levels = c("UnburntSoil","BurntSoil","PyOM"))->description1$sample
# levels(description1.melt$sample)<-c("Unburnt\nSoil","Burnt\nSoil","PyOM")
# data.table(description1.melt)->description1.melt
# description1.melt[variable=="alphaNTI",y_min := 2]  #just to set the limits of facet
# description1.melt[variable=="alphaNTI",y_max := 7]
# description1.melt[variable=="Shannon.Wiener",y_min := 5]
# description1.melt[variable=="Shannon.Wiener",y_max := 6.5]
# levels(description1.melt$variable)<-c(expression("Shannon"*"-"*"Wiener"),expression(alpha*"NTI"),"Faith",expression("Observed"~"OTUs"))
# ggplot(data=description1.melt, aes(x=sample,y=value))+
#   #  geom_segment(x=1,y=100,xend=2,yend=100,size=0.5,inherit.aes = F)+geom_segment(x=1,y=150,xend=3,yend=150,size=0.5,inherit.aes = F)+
#   geom_violin(aes(fill=sample),linetype=0)+geom_boxplot(width=0.2,lwd=0.2,outlier.size = 0.7)+
#   geom_blank(aes(y = y_min))+geom_blank(aes(y = y_max))+  #just to set the limits of facet
#   facet_wrap(~ variable, scales="free", nrow=1, labeller = label_parsed)+
#   theme_bw(base_size=10)+ guides(fill=FALSE)+
#   theme(panel.grid = element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(), 
#         axis.text.x = element_text(color="black"),
#         axis.text.y = element_text(color="black"))+
#   scale_fill_manual(values=c(color3,color1,color2))->pic1

# data.frame(sample=SampleType,pH=pH,Shannon.Wiener=shannon,alphaNTI=aNTI$NTI)->description1
# melt(description1,id=1)->description1.melt
# factor(description1$sample,levels = c("UnburntSoil","BurntSoil","PyOM"))->description1$sample
# levels(description1.melt$sample)<-c("Unburnt\nSoil","Burnt\nSoil","PyOM")
# data.table(description1.melt)->description1.melt
# description1.melt[variable=="alphaNTI",y_max := 13.5]
# description1.melt[variable=="Shannon.Wiener",y_max := 6.2]
# description1.melt[variable=="pH",y_max := 5.7]
# levels(description1.melt$variable)<-c("pH",expression("Shannon"*"-"*"Wiener"),expression(alpha*"NTI"))
# ggplot(data=description1.melt, aes(x=sample,y=value))+
#   geom_violin(aes(fill=sample),linetype=0)+geom_boxplot(width=0.2,lwd=0.2,outlier.size = 0.7)+
#   geom_blank(aes(y = y_max))+
#   facet_wrap(~ variable, scales="free", nrow=1, labeller = label_parsed)+
#   theme_bw(base_size=10)+ guides(fill=FALSE)+
#   theme(panel.grid = element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.x = element_text(color="black"),
#         axis.text.y = element_text(color="black"))+
#   scale_fill_manual(values=c(color3,color1,color2))->pic1
# 
# data.frame(sample=SampleType,pH=pH,Water=watercontent,
#            SoilOrganicCarbon=TOC,DissolvableOrganicCarbon=DOC)->description2
# melt(description2,id=1)->description2.melt
# factor(description2$sample,levels = c("UnburntSoil","BurntSoil","PyOM"))->description2$sample
# levels(description2.melt$sample)<-c("Unburnt\nSoil","Burnt\nSoil","PyOM")
# data.table(description2.melt)->description2.melt
# description2.melt[variable=="pH",y_max := 5.7]  #just to set the limits of facet
# description2.melt[variable=="Water",y_max := 0.75]
# description2.melt[variable=="SoilOrganicCarbon",y_max := 80]
# description2.melt[variable=="DissolvableOrganicCarbon",y_max := 550]
# levels(description2.melt$variable)<-c("pH","Water Content (%FW)","TOC (gC/kgDw)","DOC (mgC/kgDw)")
# ggplot(data=description2.melt, aes(x=sample,y=value))+
#   #  geom_segment(x=1,y=100,xend=2,yend=100,size=0.5,inherit.aes = F)+geom_segment(x=1,y=150,xend=3,yend=150,size=0.5,inherit.aes = F)+
#   geom_violin(aes(fill=sample),linetype=0)+geom_boxplot(width=0.2,lwd=0.2,outlier.size = 0.7)+
#   geom_blank(aes(y = y_max))+  #just to set the limits of facet
#   facet_wrap(~ variable, scales="free", nrow=1)+
#   theme_bw(base_size=10)+ guides(fill=FALSE)+
#   theme(panel.grid = element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.x = element_text(color="black"),
#         axis.text.y = element_text(color="black"))+
#   scale_fill_manual(values=c(color3,color1,color2))->pic2
# 
# data.frame(sample=SampleType,Faith=faith,Pielou=pielou,ObservedOTU=ObservedOTU,
#            SoilOrganicCarbon=TOC,DissolvableOrganicCarbon=DOC)->descriptionS1
# melt(descriptionS1,id=1)->descriptionS1.melt
# factor(descriptionS1.melt$sample,levels = c("UnburntSoil","BurntSoil","PyOM"))->descriptionS1.melt$sample
# levels(descriptionS1.melt$sample)<-c("Unburnt\nSoil","Burnt\nSoil","PyOM")
# data.table(descriptionS1.melt)->descriptionS1.melt
# descriptionS1.melt[variable=="Faith",y_max := 135]
# descriptionS1.melt[variable=="SoilOrganicCarbon",y_max := 80]
# descriptionS1.melt[variable=="DissolvableOrganicCarbon",y_max := 550]
# descriptionS1.melt[variable=="Pielou",y_max := 0.85]
# descriptionS1.melt[variable=="ObservedOTU",y_max := 1500]
# levels(descriptionS1.melt$variable)<-c(expression("Faith's"~"Diversity"),expression("Pielou's"~"J"),"Observed~OTUs",
#                                        "TOC~(gC/kgDw)","DOC~(mgC/kgDw)")
# ggplot(data=descriptionS1.melt, aes(x=sample,y=value))+
#   geom_violin(aes(fill=sample),linetype=0)+geom_boxplot(width=0.2,lwd=0.2)+
#   geom_blank(aes(y = y_max))+  #just to set the limits of facet
#   facet_wrap(~ variable, scales="free", nrow=1, labeller = label_parsed)+
#   theme_bw(base_size=12)+ guides(fill=FALSE)+
#   theme(panel.grid = element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.x = element_text(color="black"),
#         axis.text.y = element_text(color="black"))+
#   scale_fill_manual(values=c(color3,color1,color2))->picS1

data.frame(sample=SampleType,pH=pH,Shannon.Wiener=shannon,Faith=faith,Pielou=pielou,
           ObservedOTU=ObservedOTU,alphaNTI=aNTI$NTI)->description1
melt(description1,id=1)->description1.melt
factor(description1$sample,levels = c("UnburntSoil","BurntSoil","PyOM"))->description1$sample
levels(description1.melt$sample)<-c("Unburnt\nSoil","Burnt\nSoil","PyOM")
data.table(description1.melt)->description1.melt
description1.melt[variable=="alphaNTI",y_max := 13.5]
description1.melt[variable=="Shannon.Wiener",y_max := 6.2]
description1.melt[variable=="pH",y_max := 5.7]
description1.melt[variable=="Faith",y_max := 135]
description1.melt[variable=="Pielou",y_max := 0.85]
description1.melt[variable=="ObservedOTU",y_max := 1500]
levels(description1.melt$variable)<-c("pH",expression("Shannon"*"-"*"Wiener"),expression("Faith's"~"Diversity"),
                                      expression("Pielou's"~"J"),"Observed~OTUs",expression(alpha*"NTI"))
ggplot(data=description1.melt, aes(x=sample,y=value))+
  geom_violin(aes(fill=sample),linetype=0)+geom_boxplot(width=0.2,lwd=0.2,outlier.size = 0.7)+
  geom_blank(aes(y = y_max))+
  facet_wrap(~ variable, scales="free", nrow=1, labeller = label_parsed)+
  theme_bw(base_size=10)+ guides(fill=FALSE)+geom_blank(aes(y = y_max))+
  theme(panel.grid = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  scale_fill_manual(values=c(color3,color1,color2))->pic1

# pairwise.t.test(pH,SampleType,"fdr",pool.sd = FALSE)
# pairwise.t.test(shannon,SampleType,"fdr",pool.sd = FALSE)
# pairwise.t.test(faith,SampleType,"fdr",pool.sd = FALSE)
# nortest::lillie.test(pielou[isUnburnt])
# dunn.test::dunn.test(pielou,SampleType,"bh")
# pairwise.t.test(ObservedOTU,SampleType,"fdr",pool.sd = FALSE)
# pairwise.t.test(aNTI$NTI,SampleType,"fdr",pool.sd = FALSE)

dunn.test::dunn.test(pH,SampleType,"bh")
dunn.test::dunn.test(shannon,SampleType,"bh")
dunn.test::dunn.test(faith,SampleType,"bh")
dunn.test::dunn.test(pielou,SampleType,"bh")
dunn.test::dunn.test(ObservedOTU,SampleType,"bh")
pairwise.t.test(aNTI$NTI,SampleType,"fdr",pool.sd = FALSE)

wilcox.test(x=TOC[isBurnt],y=TOC[isUnburnt])
t.test(x=DOC[isBurnt],y=DOC[isUnburnt])

#betaNTI
drop.tip(root.phylotre,root.phylotre$tip.label[!(root.phylotre$tip.label %in% colnames(rare))])->root.phylotre.d
cophenetic(root.phylotre.d)->root.phylotre.d.dist
ses.beta.mntd(rare,root.phylotre.d.dist,abundance.weighted = T)->betaNTI

median(abs(betaNTI[isPyOM,isPyOM]),na.rm = T)  #median=4.89
median(abs(betaNTI[isBurnt,isBurnt]),na.rm = T)  #median=6.14
median(abs(betaNTI[isUnburnt,isUnburnt]),na.rm = T)  #median=6.71

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
  scale_fill_brewer(palette = 6,type = "qual")+ylab(expression(beta*"NTI"))+ylim(0,10)

#barplot
anno.rare[!duplicated(anno.rare[,c(2,3)]),c(2,3)]->AllPhylum
PhylumAbundance<-matrix(data=0,ncol=nrow(AllPhylum),nrow=nrow(rare))
colnames(PhylumAbundance)<-as.character(AllPhylum$Phylum)
rownames(PhylumAbundance)<-rownames(rare)
for(i in 1:nrow(AllPhylum)){
  as.matrix(rare[,as.character(anno.rare[anno.rare$Domain==AllPhylum[i,1] & anno.rare$Phylum==AllPhylum[i,2],]$OTUID)])->temp
  if(ncol(temp)>1)  rowSums(temp)->PhylumAbundance[,i] else
    temp->PhylumAbundance[,i]
}

apply(PhylumAbundance/(rowSums(PhylumAbundance)[1])>=0.01,MARGIN=2,prod)+
  apply(PhylumAbundance/(rowSums(PhylumAbundance)[1])>=0.05,MARGIN=2,sum)->temp #abundance >1% in one sample or >0.1% in all samples
as.data.frame(PhylumAbundance[,temp>=1])->PhylumAbundance.dominant
rowSums(PhylumAbundance[,temp==0])->Other
sample<-rownames(PhylumAbundance.dominant)
cbind(sample,SampleType,PhylumAbundance.dominant,Other)->PhylumAbundance.dominant

PhylumAbundance.dominant.m<-reshape2::melt(PhylumAbundance.dominant,id=c(1,2))
PhylumAbundance.dominant.m$value<-PhylumAbundance.dominant.m$value/rowSums(rare)[1]*100
PhylumAbundance.dominant.m$variable<-factor(
  PhylumAbundance.dominant.m$variable,levels=
    c("Thaumarchaeota","Acidobacteria","Actinobacteria","Bacteroidetes",
      "Chloroflexi","Cyanobacteria","Chlorophyta","V7","Planctomycetes",
      "Proteobacteria","Streptophyta","Other"))
#PhylumAbundance.bac.dominant.m$SampleType<-
#  factor(PhylumAbundance.bac.dominant.m$SampleType,levels = c("Ck","Ace","For","Lac","Pyr"))
ggplot(PhylumAbundance.dominant.m,aes(sample,value))+
  geom_col(aes(fill=variable),position = position_stack(reverse = TRUE))+
  facet_grid(~ SampleType, scales = "free_x", space = "free_x", switch = "y")+
  scale_fill_brewer(type="qual",palette = 3)+
  ylab("Abundance (%)")+
  theme_bw(12)+ 
  theme(panel.grid = element_blank(), axis.title.x=element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"))+
  guides(fill=guide_legend(title="Phylum"))->PhylumBar

#NMDS
GUniFrac(rare,root.phylotre)->unifrac
unifrac$unifracs[,,"d_0.5"]->w0.5unifrac
metaMDS(comm = w0.5unifrac,k =3)->NMDS0.5uf #stress=0.04
as.data.frame(NMDS0.5uf$points)->NMDS0.5uf.point
myfunc3(NMDS0.5uf.point,SampleType,color=c(color3,color1,color2))->ordinationAll
adonis(w0.5unifrac~SampleType,permutations = 999999) #F=11.23, p<e-6,based on factor model
adonis(w0.5unifrac~as.numeric(SampleType),permutations = 999999) #F=15.22, p<e-6,based on linear model
median(NMDS0.5uf.point$MDS1[isUnburnt])
median(NMDS0.5uf.point$MDS1[isBurnt])
median(NMDS0.5uf.point$MDS1[isPyOM])

#Shannon-NMDS1

summary(lm(faith ~ poly(NMDS0.5uf.point$MDS1,degree=2)))
summary(lm(pielou ~ poly(NMDS0.5uf.point$MDS1,degree=2)))
summary(lm(ObservedOTU ~ poly(NMDS0.5uf.point$MDS1,degree=2)))

cor.test(NMDS0.5uf.point$MDS1,pH) #r=0.77,p=4e-8
temp<-data.frame(sample=SampleType,pH=pH,NMDS1=NMDS0.5uf.point$MDS1)
ggplot(data=temp,aes(x=pH,y=NMDS1))+
  geom_point(aes(color=sample))+
  stat_smooth(method="lm", formula=y~poly(x,1), size=0.75, color="gray15")+
  theme_bw(12)+
  scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  theme(panel.grid = element_blank(), 
        legend.position = c(0.5,0.1),legend.direction= "horizontal",
        axis.text.y = element_text(color="black"))->pH.vs.NMDS1

summary(lm(shannon ~ poly(pH,degree=2)))
temp<-data.frame(sample=SampleType,Shannon=shannon,pH=pH)
ggplot(data=temp,aes(x=pH,y=Shannon))+
  geom_point(aes(color=sample))+
  geom_smooth(method="lm", formula=y~poly(x,2), size=0.75, color="gray15", linetype="dashed")+
  theme_bw(12)+  ylab("Shannon-Wiener")+
  scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  theme(panel.grid = element_blank(), 
        legend.position = c(0.5,0.1),legend.direction= "horizontal",
        axis.text.y = element_text(color="black"))->Shannon.vs.pH

summary(lm(shannon ~ poly(NMDS0.5uf.point$MDS1,degree=2)))
temp<-data.frame(sample=SampleType,Shannon=shannon,NMDS1=NMDS0.5uf.point$MDS1)
ggplot(data=temp,aes(x=NMDS1,y=Shannon))+
  geom_point(aes(color=sample))+
  stat_smooth(method="lm", formula=y~poly(x,2), size=0.75, color="gray15")+
  theme_bw(12)+  ylab("Shannon-Wiener")+
  scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  theme(panel.grid = element_blank(), 
        legend.position = c(0.5,0.1),legend.direction= "horizontal",
        axis.text.y = element_text(color="black"))->Shannon.vs.NMDS1

summary(lm(pielou ~ poly(NMDS0.5uf.point$MDS1,degree=2)))
temp<-data.frame(sample=SampleType,Pielou=pielou,NMDS1=NMDS0.5uf.point$MDS1)
ggplot(data=temp,aes(x=NMDS1,y=Pielou))+
  geom_point(aes(color=sample))+
  stat_smooth(method="lm", formula=y~poly(x,2), size=0.75, color="gray15")+
  theme_bw(12)+ylab("Pielou's J")+
  scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  theme(panel.grid = element_blank(), 
        legend.position = c(0.5,0.1),legend.direction= "horizontal",
        axis.text.y = element_text(color="black"))->Pielou.vs.NMDS1

summary(lm(ObservedOTU ~ poly(NMDS0.5uf.point$MDS1,degree=2)))
temp<-data.frame(sample=SampleType,ObservedOTU=ObservedOTU,NMDS1=NMDS0.5uf.point$MDS1)
ggplot(data=temp,aes(x=NMDS1,y=ObservedOTU))+
  geom_point(aes(color=sample))+
  stat_smooth(method="lm", formula=y~poly(x,2), size=0.75, color="gray15")+
  theme_bw(12)+ylab("Observed OTUs")+
  scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  theme(panel.grid = element_blank(), 
        legend.position = c(0.5,0.1),legend.direction= "horizontal",
        axis.text.y = element_text(color="black"))->ObservedOTU.vs.NMDS1

summary(lm(faith ~ poly(NMDS0.5uf.point$MDS1,degree=2)))
temp<-data.frame(sample=SampleType,Faith=faith,NMDS1=NMDS0.5uf.point$MDS1)
ggplot(data=temp,aes(x=NMDS1,y=faith))+
  geom_point(aes(color=sample))+
  stat_smooth(method="lm", formula=y~poly(x,2), size=0.75, color="gray15")+
  theme_bw(12)+ylab("Faith's Diversity")+
  scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  theme(panel.grid = element_blank(), 
        legend.position = c(0.5,0.1),legend.direction= "horizontal",
        axis.text.y = element_text(color="black"))->Faith.vs.NMDS1

grid_arrange_shared_legend(pH.vs.NMDS1,Shannon.vs.pH,Shannon.vs.NMDS1,Faith.vs.NMDS1,
                           Pielou.vs.NMDS1,ObservedOTU.vs.NMDS1,
                           position = "bottom",ncol = 2, nrow = 3)

#DAA
meta<-data.frame(SampleType=as.numeric(SampleType))
permute_differential_analysis(meta,t(rare),'SampleType',perm.no = 9999)->DAA.result
DAA.FDR<-data.frame(raw=DAA.result$p.raw,fdr=DAA.result$p.adj.fdr,fwer=DAA.result$p.adj.fwer)
rownames(DAA.FDR)<-names(DAA.result$p.raw)
DAA.FDR[!is.na(DAA.FDR$raw),]->DAA.FDR
dim(DAA.FDR[DAA.FDR$fdr<0.01,]->b)
apply(rare,MARGIN = 2,cor,meta,method="spearman")->rho1
apply(rare,MARGIN = 2,cor,rep(c(3,2,1),each=12),method="spearman")->rho2
which(abs(rho1)-abs(rho2)>0.01 & ifsignificant)->Sig.test
setdiff(1:2275,Sig.test)->Nul.test

#Selected OTUs
if(sum(colnames(rare)!=rownames(DAA.FDR))+sum(anno.rare$OTUID!=rownames(DAA.FDR)))
  warning("Rare, DAA.FDR, anno.rare don't match.")
which(anno.rare$OTUID %in% rownames(b))->Sig
which(!(anno.rare$OTUID %in% rownames(b)))->Nul
sign(apply(rare,MARGIN = 2,cor,meta,method="spearman"))->response
ifsignificant<-(anno.rare$OTUID %in% rownames(b))
response[!ifsignificant]<-0
if(sum(names(response)!=anno.rare$OTUID)==0)
  cbind(anno.rare,response,ifsignificant)->anno.rare else
    warning("anno.rare, response, ifsignificant don't match.")
write.csv(anno.rare[anno.rare$ifsignificant,],"SelectedOTUs.csv")

#description of selected OTUs
anno.rare$response<-as.factor(anno.rare$response)
levels(anno.rare$response)<-c("Negative","Neutral","Positive") #Positive to PyOM
length(which(anno.rare$ifsignificant)) # 853 OTUs
length(which(anno.rare$ifsignificant))/(dim(anno.rare)[1]) #37.5% of total OTUs
sum(rare[,as.character(anno.rare$OTUID[anno.rare$ifsignificant])])/sum(rare) #53.97% of total abundance
length(which(anno.rare$response=="Positive")) #525 OTUs
length(which(anno.rare$response=="Negative")) #328 OTUs

#sensitivity analysis
OTU.permutate(rare,meta$SampleType,root.phylotre,part=Sig,time=300)
OTU.permutate(rare,meta$SampleType,root.phylotre,part=Nul,time=300)
RDP<-list()
for(i in 1:(ncol(rare)%/%100)) {
  OTU.permutate(rare,meta$SampleType,root.phylotre,num=i*100,time=300)->RDP[[length(RDP)+1]]
  cat(i)
}
temp<-c()
for(i in 1:(ncol(rare)%/%100)) temp<-c(temp,RDP[[i]][2,])
temp<-matrix(temp,ncol=6,byrow = T)
temp<-as.data.frame(temp)
names(temp)<-c("Perm","UF","GUF","Euclidean","Jaccard","BC")
write.csv(temp,"RDP.csv")
temp->RDP #82.11% of importance for Sig, 20.36% of importance for Nul 

rare.ori<-rare;rare.nul<-rare;rare.sig<-rare;rare.ran<-rare
rare.nul[,Nul]<-rare[sample(1:dim(rare)[1]),Nul]
rare.sig[,Sig]<-rare[sample(1:dim(rare)[1]),Sig]
part<-sample(1:dim(rare)[2])[1:length(Sig)]
rare.ran[,part]<-rare[sample(1:dim(rare)[1]),part]
GUniFrac(rbind(rare.ori,rare.nul,rare.sig,rare.ran),root.phylotre,alpha=0.5)$unifracs[,,"d_0.5"]->GUFX4
as.data.frame(metaMDS(GUFX4,k=3,trace=F)$points)->pointX4 #stress=0.11

rare.ori.10<-rare.ori;rare.nul.10<-rare.nul
rare.sig.10<-rare.sig;rare.ran.10<-rare.ran
for(i in 1:10){
  rare.ori<-rare;rare.nul<-rare;
  rare.sig<-rare;rare.ran<-rare
  rare.nul[,Nul]<-rare[sample(1:dim(rare)[1]),Nul]
  rare.sig[,Sig]<-rare[sample(1:dim(rare)[1]),Sig]
  part<-sample(1:dim(rare)[2])[1:length(Sig)]
  rare.ran[,part]<-rare[sample(1:dim(rare)[1]),part]
  rbind(rare.ori.10,rare.ori)->rare.ori.10
  rbind(rare.nul.10,rare.nul)->rare.nul.10
  rbind(rare.sig.10,rare.sig)->rare.sig.10
  rbind(rare.ran.10,rare.ran)->rare.ran.10
}
GUniFrac(rbind(rare.ori.10,rare.nul.10,rare.sig.10,rare.ran.10),
         root.phylotre,alpha=0.5)$unifracs[,,"d_0.5"]->GUFX40
as.data.frame(metaMDS(GUFX40,k=3,trace=T,maxit=1000)$points)->pointX40

myfunc<-function(data,Sample.Type,i){
  laymat=matrix(c(rep(c(rep.int(1,3),2),times=3),rep.int(3,3),4),nrow=4,byrow = TRUE)
  if(i==1) scaleFUN <- function(x) sprintf("%.1f", x) else
    scaleFUN <- function(x) sprintf("%.2f", x)
  ggplot(data,aes(MDS1,MDS2,colour=Sample.Type,fill=Sample.Type))+
    geom_point(size=1)+
    stat_bag(prop = 1,alpha = 0.4,show.legend = T)+
    theme_bw() + theme(panel.grid = element_blank())+
    scale_color_manual(values = c(color3,"#EE6A50", "#00CD66"))+
    scale_fill_manual(values = c(color3,"#EE6A50", "#00CD66"))+
    labs(x = "NMDS1", y = "NMDS2")+xlim(-0.3,0.45)+ylim(-0.35,0.35)->pica
  ggplot(data,aes(MDS2,colour=Sample.Type,fill=Sample.Type))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw() + theme(panel.grid = element_blank(),
                       legend.position="none",
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())+coord_flip()+
    scale_color_manual(values = c(color3,"#EE6A50", "#00CD66"))+
    scale_fill_manual(values = c(color3,"#EE6A50", "#00CD66"))+
    xlim(-0.35,0.35)+ylab("Density")->picb
  ggplot(data,aes(MDS1,colour=Sample.Type,fill=Sample.Type))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw() + theme(panel.grid = element_blank(),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       legend.position="none")+scale_y_continuous(labels=scaleFUN)+
    scale_color_manual(values = c(color3,"#EE6A50", "#00CD66"))+
    scale_fill_manual(values = c(color3,"#EE6A50", "#00CD66"))+
    xlim(-0.3,0.45)+ylab("Density")->picc
  g<-ggplotGrob(pica)$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  pica<-pica+theme(legend.position="none")
  grid.arrange(pica,picb,picc,legend,layout_matrix=laymat)
}

picnes<-list()
for(i in 1:4) myfunc(pointX40[((i-1)*dim(rare)[1]*11+1):((i*dim(rare)[1])*11),],rep(SampleType,times=11),i)->picnes[[i]]
grid.arrange(picnes[[1]],picnes[[4]],picnes[[3]],picnes[[2]],ncol=2)#arranged at order "ori,ran,sig,nul"

#coherent units
anno.data<-anno.rare
response<-anno.rare$response
levels(response)<-c("Negative","NA","Positive")
anno.data$response<-response
anno.data<-anno.data[,-ncol(anno.data)]
response<-as.data.frame(response)
rownames(response)<-rownames(anno.data)
# SigClade3(anno.data[,1:7],which(anno.data$response=="Positive"),Bonferroni=0.05)->PyOM.coherent
# SigClade3(anno.data[,1:7],which(anno.data$response=="Negative"),Bonferroni=0.05)->Soil.coherent
SigClade4(anno.data[,1:7],which(anno.data$response=="Positive"),which(anno.data$response=="Negative"),Bonferroni=0.05)->all.coherent
temp<-matrix(0, nrow=nrow(all.coherent),ncol=3)
colnames(temp)<-c("PyOM","Neutral","Soil")
cbind(all.coherent,temp)->all.coherent
for(i in 1:nrow(all.coherent)) {
  anno.data[,c(as.character(all.coherent$Phylogeny)[i],"response")]->temp
  summary(temp[temp[,1]==as.character(all.coherent$taxa)[i],2])[c("Positive","NA","Negative")]->all.coherent[i,5:7]
}
write.csv(all.coherent,"CoherenceUnits.csv")

# rep(0,length(response$response))->temp1->temp2
# temp1[response$response=="Positive"]<- 1
# temp2[response$response=="Negative"]<- 1
# category<-data.frame(OTUID=anno.data$OTUID,
#                      positives=temp1,negatives=temp2,response=response$response)
# rm(temp1,temp2)

# P=rep(1,dim(Soil.coherent)[1])
# counts=matrix(rep(0,2*dim(Soil.coherent)[1]),ncol=2)
# for(i in 1:dim(Soil.coherent)[1]) try({
#   subtreeOTU=as.character(anno.data[anno.data[,as.character(Soil.coherent$Phylogeny)[i]]==as.character(Soil.coherent$taxa)[i],1])
#   dropping=setdiff(as.character(root.phylotre$tip.label),subtreeOTU)
#   subTree=drop.tip(root.phylotre,dropping)
#   phylo.d(category,subTree,"OTUID","negatives",permut = 9999)->temp
#   temp$Pval1->P[i]
#   temp$StatesTable->counts[i,]
# })
# P->OverDispersal
# colnames(counts)<-c("NullNum","SigNum")
# cbind(Soil.coherent,OverDispersal,counts)->Soil.coherent
# 
# P=rep(1,dim(PyOM.coherent)[1])
# counts=matrix(rep(0,2*dim(PyOM.coherent)[1]),ncol=2)
# for(i in 1:dim(PyOM.coherent)[1]) try({
#   subtreeOTU=as.character(anno.data[anno.data[,as.character(PyOM.coherent$Phylogeny)[i]]==as.character(PyOM.coherent$taxa)[i],1])
#   dropping=setdiff(as.character(root.phylotre$tip.label),subtreeOTU)
#   subTree=drop.tip(root.phylotre,dropping)
#   phylo.d(category,subTree,"OTUID","positives",permut = 9999)->temp
#   temp$Pval1->P[i]
#   temp$StatesTable->counts[i,]
# })
# P->OverDispersal
# colnames(counts)<-c("NullNum","SigNum")
# cbind(PyOM.coherent,OverDispersal,counts)->PyOM.coherent
# 
# 
# 
# rbind(Soil.coherent,PyOM.coherent)->CoherenceUnits
# Type<-c(rep("Soil",dim(Soil.coherent)[1]),rep("PyOM",dim(PyOM.coherent)[1]))
# cbind(CoherenceUnits,Type)->CoherenceUnits
# CoherenceUnits[,"OverDispersal"]<-p.adjust(CoherenceUnits[,"OverDispersal"],method="bonferroni")
# write.csv(CoherenceUnits,"CoherenceUnits.csv")

#phylogenetic tree
anno.rare$response->response
anno.data<-anno.rare
levels(response)<-c("Soil",NA,"PyOM")  #change level 1 to "PyOM" 
anno.data$response<-response
anno.data<-anno.data[,-ncol(anno.data)]
response<-as.data.frame(response)
rownames(response)<-as.character(anno.rare$OTUID)
# CoherenceUnits[order(CoherenceUnits$Phylogeny),]->temp3
# temp3[temp3$qvalue<0.05 & temp3$OverDispersal>0.05,]->temp
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


AnnotateTre(root.phylotre.rare,anno.data[,c(1,3)])->annotatedtre
attributes(annotatedtre[[1]])$group->group
as.character(group)->group
unique(group)->taxa
for(i in 1:length(taxa)) if(length(group[group==taxa[i]])<22) group[group==taxa[i]]<-"Other"
as.factor(group)->attributes(annotatedtre[[1]])$group

ggtree(annotatedtre[[1]],layout="fan",aes(color=group))+theme_minimal()+
  theme(panel.grid = element_blank(),legend.position="right",
        axis.text.x = element_blank(),axis.text.y = element_blank())->p
#+xlim(0, 1.8)->p
p %<+% anno.data + 
  scale_color_manual(values = c(tol21rainbow[1:12],"#D4D4D4",
                                tol21rainbow[13:19],"#888888",
                                tol21rainbow[20:23]))->p
gheatmap(p,response,colnames = FALSE,offset = -1.7, width = 0.05)->p
p+scale_fill_manual(values=c("green4","firebrick3"),na.value = NA)

#phylo.D
anno.rare$response->response
rep(0,length(response))->temp1->temp2
temp1[response=="Positive"]<- 1
temp2[response=="Negative"]<- 1
category<-data.frame(OTUID=anno.data$OTUID,
                     positives=temp1,negatives=temp2,response=response)
rm(temp1,temp2)

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
cbind(LargeClade,negatives,P)->LargeClade  #Actually, D in LargeClade is 1-D.original
rm(negatives)

phylo.d(category,phylotre,"OTUID","negatives",permut = 999)->Drecord  #1-D=0.39
phylo.d(category,phylotre,"OTUID","positives",permut = 999)->Drecord2  #1-D=0.41

#network
asin(sqrt(rare/matrix(rep(rowSums(rare),ncol(rare)),nrow=nrow(rare),byrow=F)))->rare.transformed
filter.sparsity(rare.transformed)->rare.transformed.filtered
resid(lm(rare.transformed.filtered~SampleType))->rare.transformed.filtered.residual
cor(rare.transformed.filtered.residual,method="spearman")->rho
r2p(rho,n=36)->rho.p
mat2vec(rho.p,triangle=T)->rho.p.vect
p.adjust(rho.p.vect$value,method="BH")->adj.p
mat2vec(rho,triangle=T)->rho.vect
if(sum(rho.vect$rowname!=rho.p.vect$rowname&rho.vect$colname!=rho.p.vect$colname)==0) 
  cbind(adj.p,rho.vect)->rho.vect #|rho|>0.6868, FDR<0.001
adj.rho<-matrix(as.numeric(abs(rho)>0.6868),ncol = ncol(rho))
colnames(adj.rho)<-row.names(adj.rho)<-colnames(rho)
length(which(rho.vect$value>0.6868))#1992

minerva::mine(rare.transformed.filtered.residual,n.cores=4)$MIC->MIC
MIC.ran<-list()
for(i in 1:200){
  t(apply(rare.transformed.filtered,MARGIN=1,sample))->rare.transformed.filtered.ran
  colnames(rare.transformed.filtered.ran)<-colnames(rare.transformed.filtered)
  minerva::mine(rare.transformed.filtered.ran,n.cores=3)$MIC->MIC.ran[[length(MIC.ran)+1]]
  cat(i)
}
MIC.ran.vect<-c()
for(i in 1:length(MIC.ran)){
  c(MIC.ran.vect,mat2vec(MIC.ran[[i]],triangle = T)$value)->MIC.ran.vect
}
MIC.vect<-mat2vec(MIC,triangle = T)
rm(MIC.ran)
hist(MIC)
hist(MIC.ran.vect)
rep(1,1e3)->MIC.FDR
FP=rep(1,1e3);TP=rep(1,1e3)
for(i in 1:1e3){
  FP[i]=length(MIC.ran.vect[MIC.ran.vect>(i/1e3)])
  TP[i]=length(MIC.vect$value[MIC.vect$value>(i/1e3)])
  if(i %% 100 == 0) cat(i)
}
MIC.FDR=FP/TP/length(MIC.ran.vect)*length(MIC.vect$value)
length(MIC.FDR[MIC.FDR<0.05]) #416
MIC.FDR[1e3-416] #0.05
# MIC>584, FDR<0.05
adj.MIC<-matrix(as.numeric(MIC>0.584),ncol = ncol(MIC))
colnames(adj.MIC)<-row.names(adj.MIC)<-colnames(MIC)
length(MIC.vect$value[MIC.vect$value>0.584]) #5291 edges
MIC.FDR[1000]<-0
MIC.FDR[MIC.FDR>1]<-1
MIC.vect[MIC.vect==0]<-1/1e3
cbind(MIC.vect,FDR=MIC.FDR[round(MIC.vect$value,digits = 3)*1e3])->MIC.vect
#Edges and Nodes in the network
rho.vect$colname<-as.character(rho.vect$colname)
rho.vect$rowname<-as.character(rho.vect$rowname)
MIC.vect$colname<-as.character(MIC.vect$colname)
MIC.vect$rowname<-as.character(MIC.vect$rowname)
if(sum(rho.vect$rowname!=MIC.vect$rowname & rho.vect$colname!=MIC.vect$colname)==0) cbind(MIC.vect,rho.vect)->edges
dim(edges[edges[,5]<0.001&edges[,4]<0.05,])[1]
edges[edges[,5]<0.001&edges[,4]<0.05,c(2,3,1,4,6,5,6)]->edges #FDR of spearman<0.001, that of MIC<0.05
names(edges)<-c("Source","Target","MIC","MIC.q","Spearman","Spearman.q","Weight")
sign(edges$Weight)->edges$Weight
edges$Weight*0.25+0.75->edges$Weight
write.csv(edges,"edges.csv")
nodes<-c(as.matrix(edges[,c("Source","Target")]))
nodes<-nodes[!duplicated(nodes)]
anno.data[anno.data$OTUID %in% nodes,]->nodes
colnames(nodes)[1]<-"ID"
write.csv(nodes,"nodes.csv")


#Procrustes
rownames(rare)
protest(NMDS0.5uf.point[isBurnt,],NMDS0.5uf.point[isPyOM,],permutations = 9999)
protest(NMDS0.5uf.point[isBurnt,2:3],NMDS0.5uf.point[isPyOM,2:3],permutations = 9999)
plot(protest(NMDS0.5uf.point[isBurnt,],NMDS0.5uf.point[isPyOM,],permutations = 9))
plot(protest(NMDS0.5uf.point[isBurnt,2:3],NMDS0.5uf.point[isPyOM,2:3],permutations = 9))
qplot(x=NMDS0.5uf.point[1:36,1],y=NMDS0.5uf.point[1:36,2],aes(color=SampleType),size=I(2))+
  geom_segment(aes(xend=NMDS0.5uf.point[13:24,1],yend=NMDS0.5uf.point[13:24,2],
                   x=NMDS0.5uf.point[1:12,1],y=NMDS0.5uf.point[1:12,2]),arrow=arrow(length=unit(0.2,"cm")))+
  geom_bag(aes(x=NMDS0.5uf.point[25:36,1],y=NMDS0.5uf.point[25:36,2]),prop=1,alpha=0.4,color=color3,fill=color3)+
  theme_bw(12) + theme(panel.grid = element_blank()) + xlab("NMDS1")+ylab("NMDS3")+
  scale_color_manual(values = c(color3,color1,color2))

#simulation
Specialist<-2*rare[sample(isPyOM,1e2,replace =T),Sig]+
  3*rare[sample(isUnburnt,1e2,replace =T),Sig] #Deterministic dispersal
Generalist<-5*rare[sample(c(isPyOM,isUnburnt),1e2,replace =T),Nul]
mock2<-matrix(data=0,nrow=1e2,ncol=dim(rare)[2])
mock2[,Sig]<-Specialist
mock2[,Nul]<-Generalist
mock2<-GUniFrac::Rarefy(mock2,24000)[[1]]
rownames(mock2)<-paste("Mock",1:1e2,sep ="")
rbind(rare,mock2)->mock2
vegan::diversity(mock2)->shannon2
vegan::specnumber(mock2)->OTUnum2
ses.mntd(mock2,cophenetic(root.phylotre),
         null.model ="taxa.labels",abundance.weighted=TRUE,runs = 999)->NTI.df2
NTI2<-data.frame(NTI=-NTI.df2$mntd.obs.z,SampaleType=SampleType.mock)
NTI2$NTI[isBurnt]<- -NTI.df$mntd.obs.z[isBurnt] #noticeably, NTI2$NTI[isBurnt] is identical to the one used in Fig 3
ks.test(shannon2[isMock],shannon2[isBurnt])  #D=0.46, p=0.013
ks.test(OTUnum2[isMock],OTUnum2[isBurnt])    #D=0.48, p=0.015
ks.test(NTI2$NTI[isMock],NTI2$NTI[isBurnt])  #D=0.39, p=0.056
comparison<-data.frame(SampleType=SampleType.mock,
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
  ylab("Shannon-Wiener")+ylim(5.3,6.3)+
  scale_fill_manual(values=c(color1,"gray50"))->com1
ggplot(data=comparison, aes(x=SampleType,y=OTUNum))+
    geom_violin(aes(fill=SampleType),linetype=0)+geom_boxplot(width=0.2,lwd=0.2)+
    theme_bw()+ guides(fill=FALSE)+
    theme(panel.grid = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"))+
    ylab("Observed OTUs")+ylim(1100,1500)+
    scale_fill_manual(values=c(color1,"gray50"))->com2
ggplot(data=comparison, aes(x=SampleType,y=aNTI))+
    geom_violin(aes(fill=SampleType),linetype=0)+geom_boxplot(width=0.2,lwd=0.2)+
    theme_bw()+ guides(fill=FALSE)+
    theme(panel.grid = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"))+
    ylab(expression(alpha*"NTI"))+scale_y_continuous(labels=scaleFUN,limits=c(2.2,7.5))+
    scale_fill_manual(values=c(color1,"gray50"))->com3

GUniFrac(mock2,root.phylotre)$unifracs[,,"d_0.5"]->mock2.GU
as.data.frame(metaMDS(mock2.GU,k=3)$points)->mock2.NMDS #stress=0.04
ks.test(mock2.NMDS$MDS1[isMock],mock2.NMDS$MDS1[isBurnt])  #D=0.25, p=0.43
ks.test(mock2.NMDS$MDS2[isMock],mock2.NMDS$MDS2[isBurnt])  #D=0.97, p=2e-13
ks.test(mock2.NMDS$MDS3[isMock],mock2.NMDS$MDS3[isBurnt])  #D=0.28, p=0.33
ks.test(comparison[1:12,2],comparison[13:112,2])  #shannon,D=0.46, p=0.02
ks.test(comparison[1:12,3],comparison[13:112,3])  #OTUnum,D=0.48, p=0.02
ks.test(comparison[1:12,4],comparison[13:112,4])  #OTUnum,D=0.39, n.s.


myfunc12(mock2.NMDS,SampleType.mock, c(color3,color1,"gray50",color2),-0.3,0.45,-0.35,0.25)->mockpic2
myfunc13(mock2.NMDS,SampleType.mock, c(color3,color1,"gray50",color2),-0.3,0.45,-0.25,0.25)->mockpic3
ggplot(mock2.NMDS, aes(y=MDS2, x=MDS1, colour = SampleType.mock)) +
  geom_point(size = 2,aes(alpha=SampleType.mock)) +
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_manual(name="SampleType",values=c(color3,color1,"gray50",color2)) +
  scale_alpha_manual(name="SampleType",values=c(0.3,1,1,0.3))+
  labs(y = "NMDS2", x = "NMDS1")->mockpic.legend

plots<-list(com1,com2,com3,mockpic2,mockpic3)
g <- ggplotGrob(mockpic.legend + theme(legend.position = "bottom"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$width)
layout<-matrix(c(1,2,3,4,4,5,5),nrow=1,byrow = T)
MockPic<-grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],
                      layout_matrix=layout,nrow=1)
combined<-arrangeGrob(MockPic,legend,ncol = 1,heights = unit.c(unit(1, "npc") - lheight, lheight))
grid.newpage()
grid.draw(combined)



#Fig3
margin<-15
layoutmatrix<-cbind(matrix(rep(3,7*4*margin),nrow=7),
                    rbind(cbind(matrix(rep(NA,3*4),nrow=3),
                                matrix(rep(1,3*4*(margin-1)),nrow = 3)),
                          matrix(rep(2,4*4*margin),nrow=4)))
gridExtra::grid.arrange(pic1,Shannon.vs.NMDS1,ordinationAll, layout_matrix=layoutmatrix)

#pH v.s. Shannon
index<-which(!is.na(pH))
summary(lm(shannon ~ poly(pH,degree=2)))  
summary(lm(shannon ~ pH))  
qplot(pH,shannon,color=SampleType)+
  stat_smooth(method="lm", formula=y~poly(x,2), size=0.75, color="gray15")+
  theme_bw(12)+scale_colour_manual(name= NULL,values=c(color3,color1,color2))+
  theme(panel.grid = element_blank(), 
        axis.text.y = element_text(color="black"))+
  labs(x="pH",y="Shannon-Wiener Index")+
  annotate("text", x = 4.2, y = 6.5,
           label = "paste(italic(R) ^ 2, \" = .64\")", parse = TRUE)+
  annotate("text", x = 4.2, y = 6.35,
           label = "paste(italic(P),\" = 1x\",10^-7)", parse = TRUE)->pHvsShannon
cor.test(pH,NMDS0.5uf.point$MDS1,method = "spearman") #rho=0.73, p=5e-7
qplot(pH,NMDS0.5uf.point$MDS1,color=SampleType)+
  theme_bw(12)+theme(panel.grid = element_blank(), 
                     axis.text.y = element_text(color="black"))+
  labs(x="pH",y="NMDS1")+scale_color_manual(name= "SampleType",values=c(color3,color1,color2))+
  annotate("text", y = 0.4, x = 4.16,
           label = "paste(italic(rho), \" = 0.73\")", parse = TRUE)+
  annotate("text", y = 0.35, x = 4.2,
           label = "paste(italic(P),\" = 5x\",10^-7)", parse = TRUE)->NMDS1vspH

#FigS3
layoutmatrix<-matrix(c(1,1,2,2,3,3,3,4),nrow = 2,byrow = T)
gridExtra::grid.arrange(PhylumBar,picS1,layout_matrix=layoutmatrix)
#FigS8 procrustes test for BCU-BCW, BCW-Burnt and BCU-Burnt
plot(protest1)
plot(protest2)
#FigS9 pH vs NMDS1,shannon
grid_arrange_shared_legend(pHvsShannon,NMDS1vspH,nrow=1)

qplot(x=NMDS0.5uf.point[,1],y=NMDS0.5uf.point[,2],aes(color=SampleType),size=I(2))+
  geom_segment(aes(x=NMDS0.5uf.point[1:12,1],y=NMDS0.5uf.point[1:12,2],
                   xend=NMDS0.5uf.point[13:24,1],yend=NMDS0.5uf.point[13:24,2]),arrow=arrow(length=unit(0.2,"cm")))+
  geom_bag(aes(x=NMDS0.5uf.point[25:36,1],y=NMDS0.5uf.point[25:36,2]),prop=1,alpha=0.4,color=color3,fill=color3)+
  theme_bw(12) + theme(panel.grid = element_blank(),
                       #                         legend.key.height=unit(0.35,"inch"),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank()) + 
  scale_color_manual(values = c(color3,color1,color2))

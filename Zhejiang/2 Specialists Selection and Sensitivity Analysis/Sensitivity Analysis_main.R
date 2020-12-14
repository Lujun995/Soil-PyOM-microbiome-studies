#This script should be excuted after Specialists Selection.R
source("OTU.permutate.R")
source("Post-Hoc.R")
library(GUniFrac)
library(ggplot2)
library(vegan)
library(grid)
library(gridExtra)

myfunc2<-function(data,Sample.Type){
  laymat=matrix(c(rep(c(rep.int(1,3),2),times=3),rep.int(3,3),4),nrow=4,byrow = TRUE)
  scaleFUN <- function(x) sprintf("%.1f", x)
  ggplot(data,aes(MDS1,MDS3,colour=Sample.Type,fill=Sample.Type))+
    geom_point(size=1)+
    stat_bag(prop = 1,alpha = 0.4,show.legend = T)+
    theme_bw() + theme(panel.grid = element_blank())+
    scale_color_manual(values = c("#EE6A50", "#00CD66",color3))+
    scale_fill_manual(values = c("#EE6A50", "#00CD66",color3))+
    labs(x = "NMDS1", y = "NMDS3")+xlim(-0.5,0.5)+ylim(-0.3,0.3)->pica
  ggplot(data,aes(MDS3,colour=Sample.Type,fill=Sample.Type))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw() + theme(panel.grid = element_blank(),
                       legend.position="none",
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())+coord_flip()+
    scale_fill_manual(values = c("#EE6A50", "#00CD66",color3))+
    xlim(-0.3,0.3)->picb
  ggplot(data,aes(MDS1,colour=Sample.Type,fill=Sample.Type))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw() + theme(panel.grid = element_blank(),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       legend.position="none")+scale_y_continuous(labels=scaleFUN)+
    scale_fill_manual(values = c("#EE6A50", "#00CD66",color3))+
    xlim(-0.5,0.5)->picc
  g<-ggplotGrob(pica)$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  pica<-pica+theme(legend.position="none")
  grid.arrange(pica,picb,picc,legend,layout_matrix=laymat)
}




#Calculate the pseudo-F 
#for sensitive analysis, rare.drop is used (Outliers were excluded)
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=Sig.drop,time=300)  #meta.drop$meta is numeric
OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,part=Nul.drop,time=300)
RDP<-list()
for(i in 1:(ncol(rare.drop)%/%100+1)) {
  OTU.permutate(rare.drop,meta.drop$meta,root.phylotre,num=i*100,time=300)->RDP[[length(RDP)+1]]
  cat(i)
}
temp<-c()
for(i in 1:57) temp<-c(temp,RDP[[i]][2,])
temp<-matrix(temp,ncol=6,byrow = T)
temp<-as.data.frame(temp)
names(temp)<-c("Perm","UF","GUF","Euclidean","Jaccard","BC")
temp->RDP
write.csv(temp,"RDP.csv")




#Fig S4
rare.drop.ori.10<-rare.drop.ori;rare.drop.nul.10<-rare.drop.nul;
rare.drop.sig.10<-rare.drop.sig;rare.drop.ran.10<-rare.drop.ran
for(i in 1:10){
  rare.drop.ori<-rare.drop;rare.drop.nul<-rare.drop;
  rare.drop.sig<-rare.drop;rare.drop.ran<-rare.drop
  rare.drop.nul[,Nul.drop]<-rare.drop[sample(1:dim(rare.drop)[1]),Nul.drop]
  rare.drop.sig[,Sig.drop]<-rare.drop[sample(1:dim(rare.drop)[1]),Sig.drop]
  part<-sample(1:dim(rare.drop)[2])[1:length(Sig.drop)]
  rare.drop.ran[,part]<-rare.drop[sample(1:dim(rare.drop)[1]),part]
  rbind(rare.drop.ori.10,rare.drop.ori)->rare.drop.ori.10
  rbind(rare.drop.nul.10,rare.drop.nul)->rare.drop.nul.10
  rbind(rare.drop.sig.10,rare.drop.sig)->rare.drop.sig.10
  rbind(rare.drop.ran.10,rare.drop.ran)->rare.drop.ran.10
}
GUniFrac(rbind(rare.drop.ori.10,rare.drop.nul.10,rare.drop.sig.10,rare.drop.ran.10),
         root.phylotre,alpha=0.5)$unifracs[,,"d_0.5"]->GUFX40
as.data.frame(metaMDS(GUFX40,k=3,trace=T,maxit=1000)$points)->pointX40

picnes<-list()
for(i in 1:4) myfunc2(pointX40[((i-1)*dim(rare.drop)[1]*11+1):((i*dim(rare.drop)[1])*11),],rep(Class1[c(-1,-3)],times=11))->picnes[[i]]
grid.arrange(picnes[[1]],picnes[[4]],picnes[[3]],picnes[[2]],ncol=2)#arranged at order "ori,ran,sig,nul"





library(picante)
library(ggplot2)
library(reshape2)
library(parallel)
source("ses.beta.mntd.R")

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

#betaNTI
drop.tip(root.phylotre,root.phylotre$tip.label[!(root.phylotre$tip.label %in% colnames(rare))])->root.phylotre.d
cophenetic(root.phylotre.d)->root.phylotre.d.dist
ses.beta.mntd(rare,root.phylotre.d.dist,abundance.weighted = T)->betaNTI
median(abs(betaNTI[isPyOM,isPyOM]),na.rm = T)  #median=7.45
median(abs(betaNTI[isBurnt,isBurnt]),na.rm = T)  #median=10.83
median(abs(betaNTI[isUnburnt,isUnburnt]),na.rm = T)  #median=11.81

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
  geom_violin(aes(fill=rowname),linetype=0)+geom_boxplot(width=0.2,lwd=0.2,outlier.size = 1.5)+
  theme_bw(base_size=12)+ guides(fill=FALSE)+
  theme(panel.grid = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  geom_segment(y=2,x=0.5,yend=2,xend=6.5,linetype="dashed")+
  scale_fill_brewer(palette = 6,type = "qual")+ylab(expression(beta*"NTI"))+ylim(0,17)


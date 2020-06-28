SigClade4<-function(phylotable,PyOM,Soil,Bonferroni=0.05){
  phylotablePyOM<-phylotable[PyOM,]
  phylotableSoil<-phylotable[Soil,]
  global.ratio=c(nrow(phylotablePyOM),nrow(phylotableSoil),
                 nrow(phylotable)-nrow(phylotablePyOM)-nrow(phylotableSoil))
  global.ratio=global.ratio/nrow(phylotable)
  last=ncol(phylotable)
  data=NULL
  pull_out_PyOM=rep(FALSE,nrow(phylotable))
  pull_out_Soil=rep(FALSE,nrow(phylotable))
  type="unassigned"
  
  j=0
  for(h in last:2){#count for bonferroni correction number
    taxa=levels(phylotable[,h])[unique(phylotable[,h])]
    l=length(taxa)
    for(i in 1:l){
      if(taxa[i]==""|taxa[i]=="Unassigned") next()
      if(abs(length(which(phylotablePyOM[,h]==taxa[i]))-length(which(phylotableSoil[,h]==taxa[i])))>3) j=j+1
    }
  }
  
  for(h in last:2){
    taxa=levels(phylotable[,h])[unique(phylotable[,h])]
    length(taxa)->l
    rep(1,l)->qvalue
    rep("unassigned",l)->type
    for(i in 1:l){
      if(taxa[i]==""|taxa[i]=="Unassigned") next()
      if(abs(length(which(phylotablePyOM[,h]==taxa[i]))-length(which(phylotableSoil[,h]==taxa[i])))<=3) next()
      if(sum(pull_out_PyOM[phylotable[,h]==taxa[i]])) type[i]<-"PyOM"
      if(sum(pull_out_Soil[phylotable[,h]==taxa[i]])) 
        if(type[i]=="PyOM") next() else type[i]<-"soil"  #uniqueness of type
      
      length(which(phylotablePyOM[,h]==taxa[i] & (!pull_out_PyOM[PyOM])))->n1  #sig number that is not pulled
      length(which(phylotableSoil[,h]==taxa[i] & (!pull_out_Soil[Soil])))->n2
      if((n1-n2>0)) if(type[i]=="soil") next() else type[i]<-"PyOM"
      if((n1-n2<0)) if(type[i]=="PyOM") next() else type[i]<-"soil"
      length(which(phylotable[,h]==taxa[i]))->N  #total OTU number N in the clade
      if(type[i]=="soil") length(which(phylotable[,h]==taxa[i] & (!pull_out_Soil)))->m else 
        if(type[i]=="PyOM") length(which(phylotable[,h]==taxa[i] & (!pull_out_PyOM)))->m
      if(n1+n2+m==0) if(type[i]=="PyOM") n1<-m<-1 else if(type[i]=="soil") n2<-m<-1
      #pSpecialists(N,n,p1,p2)
      pSpecialists(N,round((n1-n2)/m*N),global.ratio[1],global.ratio[2])*j->qvalue[i]
    }
    data.frame(taxa=taxa,qvalue=qvalue,type=type)->taxap
    taxap[taxap$qvalue<Bonferroni,]->sigtaxa
    Phylogeny=rep(names(phylotable)[h],time=nrow(sigtaxa))
    rbind(data,cbind(Phylogeny,sigtaxa))->data
    pull_out_Soil[which(phylotable[,h] %in% sigtaxa[sigtaxa$type=="soil",1])]<-TRUE
    pull_out_PyOM[which(phylotable[,h] %in% sigtaxa[sigtaxa$type=="PyOM",1])]<-TRUE
  }
  data
}
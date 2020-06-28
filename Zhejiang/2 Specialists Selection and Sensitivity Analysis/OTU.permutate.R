OTU.permutate<-function(otutable,group,phylotree,part="random",num=100,time=100){
  library(GUniFrac)
  library(vegan)
  library(parallel)
  n<-nrow(otutable)
  m<-ncol(otutable)
  flag<-class(part)=="character"
  
  temp=pseudo.F=pseudo.F.original=rep(0,5)
  
  #The original value
  #UniFrac
  GUniFrac(otutable,phylotree,alpha=0.5)$unifracs->unifracs
  unifracs[,,"d_0.5"]->Gunifrac
  unifracs[,,"d_UW"]->unifrac
  adonis(unifrac~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F.original[1]
  adonis(Gunifrac~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F.original[2]
  #Euclidean
  as.matrix(dist(otutable))->euclidean
  adonis(euclidean~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F.original[3]
  #Jaccard
  as.matrix(vegdist(otutable,method="jaccard"))->Jaccard
  adonis(Jaccard~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F.original[4]
  #Bray-Curtis
  as.matrix(vegdist(otutable))->BC
  adonis(BC~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F.original[5]
  
  cl <- makeCluster(detectCores()-1)  
  clusterEvalQ(cl,library(GUniFrac,vegan))
  clusterExport(cl,varlist=c("otutable","group","phylotree","part","num","n","m","flag"),envir=environment())
  parSapply(cl, 1:time, function(i,...) { 
    pseudo.F=rep(0,5)
    #randomly permutation or appointed
    if(flag) sample(1:m)[1:num]->part
    else if(class(part)=="integer") 0
    else stop("Can't identify the permutation part.")
    #permutation
    tabletemp<-otutable
    sample(1:n)->per
    tabletemp[,part]<-otutable[per,part]
    #UniFrac
    GUniFrac(tabletemp,phylotree,alpha=0.5)$unifracs->unifracs
    unifracs[,,"d_0.5"]->Gunifrac
    unifracs[,,"d_UW"]->unifrac
    adonis(unifrac~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F[1]
    adonis(Gunifrac~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F[2]
    #Euclidean
    as.matrix(dist(tabletemp))->euclidean
    adonis(euclidean~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F[3]
    #Jaccard
    as.matrix(vegdist(tabletemp,method="jaccard"))->Jaccard
    adonis(Jaccard~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F[4]
    #Bray-Curtis
    as.matrix(vegdist(tabletemp))->BC
    adonis(BC~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F[5]
    
    pseudo.F
    } )->ran
  stopCluster(cl)
  pseudo.F=apply(ran,MARGIN=1,mean)
  
  
  
  # for(i in 1:time){
  #   #randomly permutation or appointed
  #   if(flag) sample(1:m)[1:num]->part
  #   else if(class(part)=="integer") 0
  #   else stop("Can't identify the permutation part.")
  #   #permutation
  #   tabletemp<-otutable
  #   sample(1:n)->per
  #   tabletemp[,part]<-otutable[per,part]
  #   #UniFrac
  #   GUniFrac(tabletemp,phylotree,alpha=0.5)$unifracs->unifracs
  #   unifracs[,,"d_0.5"]->Gunifrac
  #   unifracs[,,"d_UW"]->unifrac
  #   adonis(unifrac~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F[1]
  #   adonis(Gunifrac~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F[2]
  #   #Euclidean
  #   as.matrix(dist(tabletemp))->euclidean
  #   adonis(euclidean~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F[3]
  #   #Jaccard
  #   as.matrix(vegdist(tabletemp,method="jaccard"))->Jaccard
  #   adonis(Jaccard~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F[4]
  #   #Bray-Curtis
  #   as.matrix(vegdist(tabletemp))->BC
  #   adonis(BC~group,permutations=1)$aov.tab$F.Model[1]->pseudo.F[5]
  #   
  #   temp=temp+pseudo.F
  # }
  if(!flag) num=length(part)
  matrix(c(0,pseudo.F.original,num,pseudo.F),nrow = 2,byrow = T)->result
  colnames(result)<-c("Perm","UF","GUF","Euclidean","Jaccard","BC")
  result
}

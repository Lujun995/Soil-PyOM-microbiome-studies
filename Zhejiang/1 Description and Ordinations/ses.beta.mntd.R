ses.beta.mntd<-function(samp, dis, abundance.weighted=FALSE, runs=99, cores="max", exclude.conspecifics=FALSE){
  if(!require("picante")) 
    stop("Install packages \"ape\", \"picante\" first.") else
      library("picante")
  library("parallel")
  comdistnt(samp,dis,abundance.weighted = abundance.weighted, exclude.conspecifics=exclude.conspecifics)->obs
  if(cores=="max"){
    cl <- makeCluster(detectCores()-1)  
    clusterEvalQ(cl,library(picante))
    clusterExport(cl,c("samp","dis", "abundance.weighted", "exclude.conspecifics"),envir=environment())
    parSapply(cl, 1:runs, function(i,...) { 
      samp->temp; 
      t(apply(temp,MARGIN=1,sample))->temp
      colnames(samp)->colnames(temp)
      comdistnt(temp,dis,abundance.weighted = abundance.weighted,
                exclude.conspecifics=exclude.conspecifics) } )->ran
    stopCluster(cl)
  } else if(cores==1){
    replicate(runs,comdistnt(temp,dis,abundance.weighted = abundance.weighted,
                             exclude.conspecifics=exclude.conspecifics))->ran
  } else if(class(cores)=="numeric"){
    cl <- makeCluster(min(cores,detectCores()-1))  
    clusterEvalQ(cl,library(picante))
    clusterExport(cl,c("samp","dis", "abundance.weighted", "exclude.conspecifics"),envir=environment())
    parSapply(cl, 1:runs, function(i,...) { 
      samp->temp; 
      t(apply(temp,MARGIN=1,sample))->temp
      colnames(samp)->colnames(temp)
      comdistnt(temp,dis,abundance.weighted = abundance.weighted,
                exclude.conspecifics=exclude.conspecifics) } )->ran
    stopCluster(cl)
  }
  apply(X = ran,MARGIN=1,mean)->ran.mean
  apply(X = ran,MARGIN=1,sd)->ran.sd
  betaNTI<-obs
  betaNTI[1:length(as.numeric(obs))]<- -(as.numeric(obs)-ran.mean)/ran.sd    #-betaMNTD, i.e. betaNTI
  as.matrix(betaNTI)->betaNTI
  diag(betaNTI)<-NA
  betaNTI
}
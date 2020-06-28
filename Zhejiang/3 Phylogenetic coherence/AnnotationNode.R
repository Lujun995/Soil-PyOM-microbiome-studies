AnnotationNode<-function(phylotable,AMBIGUITY=0.01,TOLERANT=0.20,
                         label=c("Node","node"),Missing=c("","Unassigned","unassigned")){
  #To label the nodes in the phylotable if missing vlaue>(1-AMBIGUITY) of tatal value, return "Unassigned";
  #if missing value<TOLERANT and there is one assigned value>(1-AMBIGUITY) of total assigned value, return
  #this assigned taxonomy; else return NA(not avaliable)
  #phylotable: a table contains taxonomy count of tree nodes, colnames are taxonomies with a label column,
  #            filled with total count of each taxa, datatype=(data.frame|matrix)
  #toletant: a max value to allow this node misassigned to another taxa, datatype=numeric and AMBIGUITY<0.5
  #TOLERANT: a max value to allow this node contains unassigned taxonomies, datatype=numeric and TOLERANT
  #          <1-AMBIGUITY
  #label: the colnames of the label column, datatype=character
  #Missing: the colnames of the unassigned value column, datatype=character
  
  #value: return phylotable with an representative taxonomy column, datatype=data.frame
  
  
  if(!(class(phylotable) %in% c("data.frame","matrix")))
    stop("Unknow datatype of \"phylotable\"")
  if(AMBIGUITY>=0.5) warning("risk of ambigously assignment because of a larger \"AMBIGUITY\"")
  if(TOLERANT>=1-AMBIGUITY) warning("risk of ambigously assignment because of a larger \"TOLERANT\"")
  
  colnames(phylotable)->taxa
  which(!(taxa%in%label|taxa%in%Missing))->valid
  which(taxa%in%Missing)->invalid
  length(valid)->n
  length(invalid)->m
  if(m>1){
    rowSums(phylotable[,invalid])->invalidsum
    cbind(phylotable[,valid],invalidsum)->data
    m=1
  }else 
    if(m==1) 
      cbind(phylotable[,valid],phylotable[,invalid])->data else 
        if(m==0) phylotable[,valid]->data
  as.matrix(data)->data
  if(m==1) colnames(data)<-c(taxa[valid],"Unassigned") else
    colnames(data)<-c(taxa[valid])
  taxa[valid]->validtaxa
  
  rowSums(data)->total
  if(m==1) invalidratio<-data[,n+m]/total
  rowSums(data[,1:n])->validtotal
  matrix(rep(validtotal,times=n),ncol=n,byrow = FALSE)->validtotal
  validratio<-data[,1:n]/validtotal
  
  represent<-rep("NA",time=nrow(phylotable))
  
  if(m==1) {
    represent[which(invalidratio>=(1-AMBIGUITY))]<-"Unassigned"
    a=(validratio>1-AMBIGUITY)
    b=(invalidratio<=TOLERANT)
    b=matrix(rep(b,times=n),ncol=n,byrow = FALSE)
    c=a&b
    for(i in 1:n) represent[which(c[,i])]<-validtaxa[i]
  }
  if(m==0) {
    a=(validratio>1-AMBIGUITY)
    for(i in 1:n) represent[which(a[,i])]<-validtaxa[i]
  }
  
  colnames(phylotable)[which(colnames(phylotable)=="")]<-"Unassigned"
  phylotable<-as.data.frame(phylotable)
  return(cbind(phylotable,represent))
}
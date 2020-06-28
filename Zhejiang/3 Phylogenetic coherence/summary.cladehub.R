summary.cladehub<-function(cladehub){
  rowSums(cladehub[,2:(ncol(cladehub)-1)])->total
  data.frame(Node=cladehub[which(total>100),1],
             Represent=cladehub[which(total>100),ncol(cladehub)],
             Number=total[which(total>100)])->temp
  return(temp)
}
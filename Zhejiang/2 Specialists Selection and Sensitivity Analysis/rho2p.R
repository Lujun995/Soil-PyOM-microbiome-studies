rho2p<-function(rho,observed_n){
  #a more robust and precise function to calculate the p-value in pearson 
  #correlation matrix comparing to psych::corr.test
  rho=abs(as.matrix(rho))
  nf=observed_n-2
  t=(nf*(1/(1-rho^2)-1))^0.5
  pt(c(t),nf,lower.tail=FALSE)*2->p
  matrix(p,nrow=nrow(rho),ncol=ncol(rho))->p
  rownames(rho)->rownames(p)
  colnames(rho)->colnames(p)
  p
}
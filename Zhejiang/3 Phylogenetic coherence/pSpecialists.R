pSpecialists<-function(N,n,p1,p2){
  #p1 is the probability of positive
  #p2 is the probability of negative
  #N is the sample size
  #n is the obeserved sum
  #two-sided possibility
  sum=0
  if(n<0) {
    p1->p
    p2->p1
    p->p2
    n<- -n
  }
  for(i in 0:((N-n)%/%2)) { #i: the number of possible negatives
    for(j in (n+i):(N-i)) {#j: the number of possible positives
      sum=sum+dmultinom(c(j,i,N-i-j),prob=c(p1,p2,1-p1-p2))
    }}
  sum*2
}
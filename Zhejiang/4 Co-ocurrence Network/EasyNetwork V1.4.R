#EasyNetwork is useful for people who want to do something with 
#Network Analysis though with limited experience. Enjoy the convience! 
#Copyright (C) 2017  Chung Loojuen
 
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
 
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
 
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see http://www.gnu.org/licenses/.

setup<-function(){
#to set up the Bioconductor
	if(!require(multtest)) {
		source("http://bioconductor.org/biocLite.R")
		options("BioC_mirror"="http://mirrors.ustc.edu.cn/bioc/")
		biocLite()
		library(multtest)
	};
}

mat2vec<-function(x,triangle=FALSE,data.in=NULL){
#to convert a matrix or data.frame to three vectors: value,rname,cname
#especially when the matrix or data.frame is symmetric or triangle.
	rownames(x)->rnam
	colnames(x)->cnam
	if(!is.matrix(x)){
		if(is.data.frame(x)) as.matrix(x)->x else 
		  stop("Input is not a dataframe or matrix.")
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
		if(data.in=="lower") x[upper.tri(x,diag=TRUE)]<-NA else 
		  if(data.in=="upper") x[lower.tri(x,diag=TRUE)]<-NA else 
        stop("You said it's a triangle matrix, but where is
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

r2p<-function(x,n){
#a more robust and precise function to calculate the p-value in pearson 
#correlation matrix comparing to psych::corr.test
	x=as.matrix(x)
	nf=n-2
	t=(nf*(1/(1-x^2)-1))^0.5
	pt(c(t),nf,lower.tail=FALSE)*2->p
	matrix(p,nrow=nrow(x),ncol=ncol(x))->p
	rownames(x)->rownames(p)
	colnames(x)->colnames(p)
	p
}

function(r,df) sqrt(df*(1/(1-r^2)-1))*sign(r)

table2p<-function(value,table){
#to convert correlation score matrix to a p-value according to existed 
#value_p table and the very firstcolumn in the table should be 
#correlation scores and the second should be the corresponding p_value.
#NOTE: Suppose: Higher score with lower p.

        table=table[order(table[,2]),]
        lt=nrow(table)
        value<-as.matrix(value)
        x=rep(1,times=nrow(value)*ncol(value))

        x[which(value>=table[1,1])]<-table[1,2]
        for(i in 2:lt){ 
#    x[which(value>=table[i,1]&value<table[i-1,1])]<-table[i,2]
        }
        x<-matrix(x,nrow(value),ncol(value))
        rownames(x)<-rownames(value)
        colnames(x)<-colnames(value)
x
}

filter.abundance<-function(otu,threshold=0.0001){
#to filter OTUs that have less than "threshold" relative abundance
  apply(otu,MARGIN = 2,sum)->x
  otu[,which(x>=(threshold*sum(x)))]
}

filter.sparsity<-function(otu,tol=0.4){
#to filter OTUs that have more than "tol"erant sparsity
  n=nrow(otu)
  sparsity=colSums(otu==0)
  otu[,which(sparsity<=(tol*n))]
}

wipe.r<-function(net,r=0.8){
  #Calculate RMT threshold in http://129.15.40.240/mena/ first.
  #The defualt value is just for test.
  net[abs(net$r)>=r,]
}

EasyNetwork<-function(otu,rvalue=0.8,alpha=0.01,method="FDR"){
#This function is designed for naive users who want to get started in Network 
#Analysis. The Spearman's correlation is used to calculate Network edges, and 
#the multitest is controlled by either FDR(BH method) or FWER(Bonferroni 
#method). OTUs that have more than 40% sparsity are filtered
#(Sophie Weiss,2016).
#Calculate RMT threshold in http://129.15.40.240/mena/ first, and pass it to 
#rvalue.
  
#to filter OTUs
#	wipe.abundance(otu)->otu #Uncomment it if rela. abun. should be filtered.
	filter.sparsity(otu)->otu
  nr=nrow(otu)
  #nc=ncol(otu)
#to calculate the p&r                                          
	stats::cor(otu,method="spearman")->r
	r2p(r,nr)->p
#to control multi-test
	mat2vec(p,tri=TRUE,data.in="lower")->p
	mat2vec(r,tri=TRUE,data.in="lower")->r
	if(identical(p[,2],r[,2])&identical(p[,3],r[,3]))
		data<-data.frame(p=p[,1],r=r[,1],source=p[,2],target=p[,3])
	else stop("I hope this massage will never show...")
	if(method=="FDR") proc="BH"
	else if(method=="FWER") proc="Bonferroni"
	else stop("Method not in (FDR|FWER).")
	multtest::mt.rawp2adjp(p$value,proc,alpha)->ADJ;
	adjp<-ADJ$adj[order(ADJ$index),]
	if(identical(adjp[,1],data$p)){
	    if(method=="FDR") {cbind(data,FDR=adjp[,2])->data; data[data$FDR<alpha, ]->data;}
	    else if(method=="FWER") {cbind(data,FWER=adjp[,2])->data; data[data$FWER<alpha, ]->data;}
	}
	else stop("adjp[,1]!=data$p")
	data=wipe.r(data,rvalue)
	rownames(data)<-1:nrow(data)
#to format 
	type<-rep("undirected",times=nrow(data))
	weight<- -log10(data$p+1e-30)
  cbind(data,type,weight)->data
#   min(abs(data$r))->r.min
#   max(abs(data$r))->r.max
#   color=rep(0,times=nrow(data))
#   cbind(data,color)->data
#   data[data$r>0,]$color<-(data[data$r>0,]$r-r.min)/(r.max-r.min)
#   data[data$r<0,]$color<-(data[data$r<0,]$r+r.min)/(r.max-r.min)
data
}


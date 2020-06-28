#incomplete!!!!!
read.csv("16S.uparse.csv",header=TRUE)->otudata
rownames(otudata)<-otudata[,1]
otudata<-otudata[,-1]
data<-otudata[,-ncol(otudata)]
data<-data[,-((ncol(data)-2):ncol(data))]
colnames(data)<-c("Unburnt1","Unburnt2","Unburnt3","Unburnt4","Unburnt5","Unburnt6",
                  "Unburnt7","Unburnt8","Unburnt9","Unburnt10","Unburnt11","Unburnt12",
                  "BCU1_2","BCU1_3","BCU2_1","BCU2_2","BCU3_1","BCU3_3","BCU4_1","BCU4_2",
                  "BCU4_3","BCW1_1","BCW1_2","BCW1_3","BCW2_1","BCW2_2","BCW2_3","BCW3_1",
                  "BCW3_2","BCW3_3","BCW4_1","BCW4_2","BCW4_3","Burnt1_1","Burnt1_2",
                  "Burnt1_3","Burnt2_1","Burnt2_2","Burnt2_3","Burnt3_1","Burnt3_2",
                  "Burnt3_3","Burnt4_1","Burnt4_2","Burnt4_3")
read.csv("rare.csv",header=TRUE)->rare
rownames(rare)<-rare[,1]
rare<-rare[,-1]
rare<-rare[,-((ncol(rare)-2):ncol(rare))] #remove the last three columns
colnames(rare)<-colnames(data)
rare.drop<-rare[,c(-1,-3)]
t(rare)->rare
t(rare.drop)->rare.drop
data<-t(data)
read.csv("Phylogeny.csv",header=TRUE)->phylo
ape::read.tree("16S.refined.tre")->phylotre
ape::root(phylotre, 1, r = TRUE)->root.phylotre

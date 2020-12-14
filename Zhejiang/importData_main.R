#Prepare the OTU table, the taxonomy table, the phylogenetic tree 
#and meta data for further analyses.
#You can directly use the .RData for convenience. 
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
Sample.Type<-as.factor(c(rep("UnburntSoil",12),rep("PyOM",21),rep("BurntSoil",12)))
PyOM.Type<-c(rep("Unwashed",9),rep("Washed",12))
#meta is the PyOM scale with 1 = unburnt soils, 2 = burnt soils and 3 = PyOM.
#However, meta here is a soil scale with 3 = UnburntSoil, 1 = PyOM, 2 = BurntSoil; therefore,
#in the following analyses, OTU response should be reversed (-1 = positive to PyOM effect).
meta <- c(rep(3,12),rep(1,21),rep(2,12)) 
meta <- data.frame(meta=meta)
Class1<-Sample.Type
Class2<-PyOM.Type
color1="#EE6A50"
color2="#00CD66"
color3="#5CACEE"
isBurnt=34:45
isPyOM=13:33
isUnburnt=1:12
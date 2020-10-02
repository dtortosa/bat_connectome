###library
library(ape)

### read the distance matrix calculated with HGC_matrix.py
m <- as.matrix(read.table("/Users/diegosalazar/Google Drive/science/open_projects/connectome/code/connectome/tree/distance_matrix.txt", head=T, row.names=1))

### aply algorithm of nearest neighbourh on matrix to create the arbol
arbol <- nj(as.dist(m))
arbol$edge.length <- arbol$edge.length + abs(min(arbol$edge.length))+1

### read data of function related with BAT of all UCP1 connectome genes
func = read.csv("/Users/diegosalazar/Google Drive/science/open_projects/connectome/code/connectome/tree/relation_with_BAT.csv", header=T)

### color of tip labels according to BAT function 
colr<-rep("black", nrow(func)) #create a vector with length equal to the number of genes, all black
colr[func$FUN == 1]<-"red" # in red all genes related with BAT
colr[func$FUN == 2]<-"blue" # in blue proteinas relacionadas con otras proteínas que se relacionana a su vez con el BAT, pero no hay descrita una relacion directa
colr[func$FUN == 0]<-"blue" #in black novel candidate genes
names(colr)<-func$GENES #genes as names of vector
colr<-colr[arbol$tip.label] #order col vector with tip label position in the tree

## plot the tree
pdf("/Users/diegosalazar/Google Drive/science/open_projects/connectome/code/connectome/tree/FGA_tree.pdf",width=10, height=10,pointsize=0.1)
plot(arbol,cex=1.1,edge.width=1, font=15, type="fan", tip.color=colr)
#legend("topright", legend=c("Direct relation described", "Indirect relation described", "No relation described"), fill= c("red", "blue", "black"), cex=2)
dev.off()

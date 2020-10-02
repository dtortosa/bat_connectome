###Figure of biological distance of UCP genes and a sample of human genes (1 million)

####load data of all the pairs across genome by menad of a RData (faster)
load("/SCRATCH/UGR002/dsalazar/ucp1_connectome/data/connectome_data_loaded.RData")

#############################################
######## READ THE UCP1 CONNECTOME ###########
#############################################
ucp1_conn = read.table("/SCRATCH/UGR002/dsalazar/ucp1_connectome/data/UCP1.txt", sep="\t", header=T)
str(ucp1_conn)
head(ucp1_conn)


#############################################
######## SELECT TOP 1% BY P.VALUE ###########
#############################################

##Select the top 1% of all human genes by p-value, only genes close to the core by a p.value lower than 0.01. You have to use the core p.value, which is the p.value about the distance between the core of your connectome and the target gen. The next column (target to core), is the p.value from the target to the core gene but in the connectome of the target gene. For example, Target_in_source_P.value.percentile. of UCP2 in the connectome of UCP1 is 0.00006, whilst the Source_in_target_P.value.percentile. is 0.00012, exactly the Target_in_source_P.value.percentile. of UCP1 in the connectome of UCP2. 
subset = ucp1_conn[ucp1_conn$Target_in_source_P.value.percentile<0.01,]
subset$Target
length(subset$Target)

################################################################
######## distance between UCP1 conectome genes #################
################################################################

#####Loop for createing a vector with all possible pairs of genes included en in the top 1% of UCP1 connectome.
ucp_conn_pairs = NULL
for(i in 1:nrow(subset)){ #for each row of subset
    selected_row = subset[i,] #select the row [i]
    if(!i==nrow(subset)){ #if the it is not the final row
        subset_2 = subset[(i+1):168,] #select all the rows following the [i] row, in this way each pair is only written one time      
        for(j in 1:nrow(subset_2)){ #for each gen of subset2
            selected_row_2 = subset_2[j,] #select the gen [j]
            ucp_conn_pairs = append(ucp_conn_pairs, paste(selected_row$Target, selected_row_2$Target, sep="_")) #create a string with the [i] gen and the [j] gen.
        }
        print(nrow(subset_2)) #print nrows to confirm that each step reduce the number of rows in 1. 
    }
}

str(ucp_conn_pairs)
head(ucp_conn_pairs)

#tests
head(ucp_conn_pairs, 168) #begins UCP2 with UCP3, not UCP1
head(ucp_conn_pairs, 334) #begins UCP3 with BMP7, not UCP2

####create a colum with pairs between genes
all_pairs$pairs = paste(all_pairs$gen1, all_pairs$gen2, sep="_")

####select those rows form all_pairs in which the pair is included in the UCP1 connectome
pairs_ucp1_connec = all_pairs[which(all_pairs$pairs %in% ucp_conn_pairs),]

###set pseudo random number
seed = 4675
set.seed(seed)

#### sample of all genes (1 million)
pairs_random_connec = all_pairs[sample(1:nrow(all_pairs), 1000000),]


#### calculate proportions of pairs in each categoiy of biological distance
proportions = rbind.data.frame( #number of pairs in each category divided by the total of pairs. We will do for UCP1 and sample genes
    nrow(pairs_ucp1_connec[pairs_ucp1_connec$biological_distance < 10,])/nrow(pairs_ucp1_connec),
    nrow(pairs_ucp1_connec[pairs_ucp1_connec$biological_distance >= 10 & pairs_ucp1_connec$biological_distance < 20,])/nrow(pairs_ucp1_connec),
    nrow(pairs_ucp1_connec[pairs_ucp1_connec$biological_distance >= 20 & pairs_ucp1_connec$biological_distance < 30,])/nrow(pairs_ucp1_connec),
    nrow(pairs_ucp1_connec[pairs_ucp1_connec$biological_distance >= 30 & pairs_ucp1_connec$biological_distance < 40,])/nrow(pairs_ucp1_connec),
    nrow(pairs_ucp1_connec[pairs_ucp1_connec$biological_distance >= 40,])/nrow(pairs_ucp1_connec),
    nrow(pairs_random_connec[pairs_random_connec$biological_distance < 10,])/nrow(pairs_random_connec),
    nrow(pairs_random_connec[pairs_random_connec$biological_distance >= 10 & pairs_random_connec$biological_distance < 20,])/nrow(pairs_random_connec),
    nrow(pairs_random_connec[pairs_random_connec$biological_distance >= 20 & pairs_random_connec$biological_distance < 30,])/nrow(pairs_random_connec),
    nrow(pairs_random_connec[pairs_random_connec$biological_distance >= 30 & pairs_random_connec$biological_distance < 40,])/nrow(pairs_random_connec),
    nrow(pairs_random_connec[pairs_random_connec$biological_distance >= 40,])/nrow(pairs_random_connec))

###bind proportions and gene and distance classes
data = cbind.data.frame(
    c(rep("UCP1 connectome genes", 5), rep("All genes", 5)),    
    rep(c("<10", "<20", "<30", "<40", ">=10"),2),
    proportions)      
colnames(data) <- c("gene_class", "distance_class", "proportions")
data

###Creamos una matriz a #partir de los datos.
xtabs.data<-xtabs(proportions ~factor(gene_class)+factor(distance_class), data=data)
xtabs.data

###plot
cairo_pdf(paste("/SCRATCH/UGR002/dsalazar/ucp1_connectome/results/barplot_proportions_", seed, ".pdf", sep=""))#we use cairo because pdf() gives problems with \u2265 (simbol >=)

#barplot
xs<-barplot(xtabs.data, beside=TRUE, xpd=F, axes=FALSE, axisnames=F, ylim=c(0,1), width=3, lwd=1:2, col=c("gray", "red"))

#axis texts
mtext(text=expression(bold("")), side = 3, line=1.8, outer=FALSE, cex=0.8)
mtext(text=expression(bold("Proportion")), side = 2, line=2.8, outer=FALSE, cex=0.8)
mtext(text=expression(bold("Biological distance")), side = 1, line=3, outer=FALSE, cex=0.8)

#axis
axis(side=1, at=c(-10, xs[1,1],xs[1,2], xs[1,3], xs[1,4], xs[1,5], 42) ,labels=c("","", "","", "", "", ""), cex.lab=0.9, cex.axis=1, lwd=1, lwd.ticks=0, font=2, mgp=c(3, 1, 0))
axis(side=2, cex.lab=1.6, cex.axis=1, font=2, mgp=c(3, 1, 0))

#axis labels of X axis
mtext(at=(xs[1,1]+xs[2,1])/2, text=expression(bold("<10")), side = 1, line=1, outer=FALSE, cex=0.8)
mtext(at=(xs[1,2]+xs[2,2])/2, text=expression(bold("<20")), side = 1, line=1, outer=FALSE, cex=0.8)
mtext(at=(xs[1,3]+xs[2,3])/2, text=expression(bold("<30")), side = 1, line=1, outer=FALSE, cex=0.8)
mtext(at=(xs[1,4]+xs[2,4])/2, text=expression(bold("<40")), side = 1, line=1, outer=FALSE, cex=0.8)
mtext(at=(xs[1,5]+xs[2,5])/2, text=expression(bold("\u2265"~"40")), side = 1, line=1, outer=FALSE, cex=0.8)

#legend
legend(x=xs[2,4], y=1, legend=c(expression(bold("All genes")), expression(bold("UCP1 genes"))), cex=0.9, xpd=T, fill=c("gray", "red"))
dev.off()

#eliminate the bid file with all pairs of genes
rm(all_pairs)

#save the results
save.image(paste("/SCRATCH/UGR002/dsalazar/ucp1_connectome/results/test3_ucp_figure_", seed, ".RData", sep=""))
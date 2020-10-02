############################
########## TABLE 1 #########
############################

### read data of function related with BAT of all UCP1 connectome genes
func = read.csv("/Users/diegosalazar/Google Drive/science/open_projects/connectome/code/connectome/tree/relation_with_BAT.csv", header=T)
func$FUN = as.factor(func$FUN)
str(func)


### read UCP1 connectome
ucp1 = read.table("/Users/diegosalazar/Google Drive/science/open_projects/connectome/data/human_connectome/UCP1.txt", sep="\t", header=T)
str(ucp1)

### select genes reported to be related with BAT
kown_genes = ucp1[which(ucp1$Target %in% func[which(func$FUN == 1),]$GENE),]
str(kown_genes)

### select columns of interest 
table_1 = kown_genes[,which(colnames(kown_genes) %in% c("Target", "Distance", "Rank", "Target_in_source_P.value.percentile.", "Route"))]

### change column names 
colnames(table_1)[which(colnames(table_1) == "Target")] <- "Genes"
colnames(table_1)[which(colnames(table_1) == "Distance")] <- "Biological distance between known and core"
colnames(table_1)[which(colnames(table_1) == "Rank")] <- "Ranks of known in core"
colnames(table_1)[which(colnames(table_1) == "Target_in_source_P.value.percentile.")] <- "p-Value (percentile) of known in core"
colnames(table_1)[which(colnames(table_1) == "Route")] <- "Route between known and core"

### drop UCP1 row
table_1  = table_1[-1,]
str(table_1)

###save table
write.table(table_1, "/Users/diegosalazar/Google Drive/science/open_projects/connectome/results/connectome_results/tables/table_1.csv", sep=",", col.names = TRUE, row.names = FALSE)

############################
####### APPENDIX S1 ########
############################
table_2 = func[order(func$GENES),]

table_2$FUN[which(table_2$FUN == 2)] <- 0

colnames(table_2)[which(colnames(table_2)=="GENES")] <- "Genes"
colnames(table_2)[which(colnames(table_2)=="FUN")] <- "BAT relationship"

###save table
write.table(table_2, "/Users/diegosalazar/Google Drive/science/open_projects/connectome/results/connectome_results/tables/appendix_S1_ordered.csv", sep=",", col.names = TRUE, row.names = FALSE)





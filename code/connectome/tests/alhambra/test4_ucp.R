#############################################
######## READ THE UCP1 CONNECTOME ###########
#############################################
ucp1_conn = read.table("/SCRATCH/UGR002/dsalazar/ucp1_connectome/data/UCP1.txt", sep="\t", header=T)
str(ucp1_conn)
head(ucp1_conn)

### create a vector with all codifican human genes ordered alphabetically
all_genes = sort(ucp1_conn$Target)

### extract top 1% of UCP1 connectome 
subset = ucp1_conn[ucp1_conn$Target_in_source_P.value.percentile<0.01,]
length(subset$Target)

### read data of function related with BAT of BAT connectome
func = read.csv("/SCRATCH/UGR002/dsalazar/ucp1_connectome/data/relation_with_BAT.csv", header=T)
func$FUN = as.factor(func$FUN)
str(func)

### drop UCP1 because UCP1 was used as core genes, thus of course will be related with BAT, we used because of this! Therefore, we will have a list of 167 genes related with according to the bibliography
func = func[-which(func$GENE == "UCP1"),]
subset = subset[-which(subset$Target == "UCP1"),]

### from the list of all genes select those which are relate with BAT according to bibliography and are included in the BAT connectome 
BAT_described = all_genes[which(all_genes %in% func[which(func$FUN==1),]$GENES)]

### loop to run the simulations
seed = 1
n_iter = 1000000
set.seed(seed)

results = NULL
for(i in 1:n_iter){

    ### random selection of as many genes as those described to be related with BAT
    random_genes = all_genes[sample(1:length(all_genes), length(BAT_described))]

    ### from random genes select those included in the BAT connectome 
    random_in_BAT_connectome = random_genes[which(random_genes %in% subset$Target)]

    ### if the number of random genes included in the BAT connecotme is the same than the list of BAT-described genes included in the connectome 
    if(length(random_in_BAT_connectome) == length(BAT_described)){
        results = append(results, "YES")
    } else {
        results = append(results, "NO")
    }  
}

### save results
write.table(results, "/SCRATCH/UGR002/dsalazar/ucp1_connectome/results/results_test_4.txt", col.names = FALSE, row.names = FALSE, sep="\t")

#Calculate p.vale as the proability of a random pair of genes will have lower biological distance than the median of biological distance between all UCP1 connectome genes.
pval = prop.test(x=length(which(results == "YES")), n=length(results)) #https://stats.stackexchange.com/questions/167164/test-of-equal-proportions-with-zero-successes

#print pval 
pval
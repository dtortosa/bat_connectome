##############################################################################################
####################################### TEST 6 ############################################### 
##############################################################################################

#Test that the Biological distance between known BAT genes and unknown BAT genes inside the connectome are smaller than that of known BAT and random genes.

#we take two BAT genes, one known to be associated with BAT and another not known to be associated with BAT. We compare their distance with the distance to a random gene. This is the comparison.



#############################################
######## LOAD DATA OF ALL THE PAIRS #########
#############################################

#load data of all the pairs across genome, i.e., the biological distance between the 16000 genes included in the conectome.
#YOU NEED TO CONNECT THE HD!!
load("/media/dftortosa/easystore/msi_diego_24_11_2019/science_big_documents/UCP_connectome/connectome_data_loaded.RData")



#############################################
######## READ THE UCP1 CONNECTOME ###########
#############################################

#load the connectome with UCP1 as core gene
ucp1_conn = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/open_projects/connectome/data/human_connectome/UCP1.txt", sep="\t", header=T)
str(ucp1_conn)
head(ucp1_conn)



#############################################
######## READ THE BAT INFORMATION ###########
#############################################

#load the information about BAT relationships
bat_relationship = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/open_projects/connectome/results/connectome_results/tables/appendix_S1_ordered.csv", sep=",", header=TRUE)
str(bat_relationship)
summary(bat_relationship)
head(bat_relationship)


###########################################
######## CREATE GROUPS OF GENES ###########
###########################################
#create samples of known and unknown BAT genes within the connectome, and random genes

#select those genes that are known to be associated with BAT within the connectome
known_bat_genes = bat_relationship[which(bat_relationship$BAT.relationship == 1),]$Genes

#select those genes that are NOT known to be associated with BAT within the connectome
unknown_bat_genes = bat_relationship[which(bat_relationship$BAT.relationship == 0),]$Genes


## select all genes outside the connectome (both known and unknown to be associated with BAT)
#create a vector with all codifican human genes ordered alphabetically from the UCP1 connectome
all_genes = sort(ucp1_conn$Target)
#select all genes except those included in the connectome
random_genes = all_genes[-which(all_genes %in% bat_relationship[which(bat_relationship$BAT.relationship %in% c(0,1)),]$Genes) ]


#check length
length(known_bat_genes)
length(unknown_bat_genes)
length(random_genes)
length(known_bat_genes) + length(unknown_bat_genes) + length(random_genes) == nrow(ucp1_conn) #we have all the genes included considered in the UCP1 connectome



#############################################
######## SAMPLE EACH GROUP AND BIND #########
#############################################

#select the sample size
sample_size_test = 100000

#set the seed to have reproducible results
set.seed(45788)

#randomly select genes from the three groups. The total number is "sample_size_test". As the the known and unknown sets are very small compare to the rest of the genome, we need replace=TRUE to have the possibility to select the same gene several times. 
sample_known_bat_genes = known_bat_genes[sample(1:length(known_bat_genes), sample_size_test, replace=TRUE)]
sample_unknown_bat_genes = unknown_bat_genes[sample(1:length(unknown_bat_genes), sample_size_test, replace=TRUE)]
sample_random_genes = random_genes[sample(1:length(random_genes), sample_size_test, replace=TRUE)] #replace let you to select the same number several times. 
    #But what if a value can be selected multiple times? This is known as sampling with replacement. Sample supports this via an additional parameter: replace. https://www.programmingr.com/examples/neat-tricks/sample-r-function/
    #If replace=FALSE, each iteration you have less sample where take data.

#bind all the data. This will be the genes to be compared.
genes_to_test = cbind.data.frame(sample_known_bat_genes, sample_unknown_bat_genes, sample_random_genes)
nrow(genes_to_test) == sample_size_test


##create a row index variable to have ID for each row
#create the variable
genes_to_test$row_index = 1:nrow(genes_to_test)
#check the variable is correct
summary(genes_to_test$row_index == 1:nrow(genes_to_test))
#reorder columns
genes_to_test = genes_to_test[,c("row_index", "sample_known_bat_genes", "sample_unknown_bat_genes", "sample_random_genes")]


##independently to the other, you got the same biological distance
all_pairs[which(all_pairs$gen1 == "UCP1" & all_pairs$gen2 == "FTO"),]
all_pairs[which(all_pairs$gen1 == "FTO" & all_pairs$gen2 == "UCP1"),] #the same biological distance
    #HABRIA QUE COGER UNA MUESTRA AL AZAR Y COMPRAR QUE DA IGUAL EL ORDEN DEL PAR, LA DISTANCIA BIOLOGICA ES LA MISMA..



###################################
######## SUBSET ALL PAIRS #########
###################################

#REVISA ESTA SECCION

#all pairs is a very big file (more than hundred million rows), so we are going to select only those rows we need

#select those pairs that have as gen 1 a known bat gene and as gene2 a unknown bat gene
all_pairs_sample_1 = which(all_pairs$gen1 %in% as.vector(genes_to_test$sample_known_bat_genes) & all_pairs$gen2 %in% as.vector(genes_to_test$sample_unknown_bat_genes))

#select those pairs that have as gen 1 a known bat gene and as gene2 a random gene
all_pairs_sample_2 = which(all_pairs$gen1 %in% as.vector(genes_to_test$sample_known_bat_genes) & all_pairs$gen2 %in% as.vector(genes_to_test$sample_random_genes))

#bind the rows
all_pairs_sample = c(all_pairs_sample_1, all_pairs_sample_2)

#remove duplicates
all_pairs_sample = all_pairs_sample[which(!duplicated(all_pairs_sample))]

#select these rows from all pairs
all_pairs_subset = all_pairs[all_pairs_sample,]

#remove all pairs
#rm(all_pairs)



###################################
######## MAKE THE TEST6 ###########
###################################

#required package
require(plyr)

#write a function to make the test 6
#test6_row = genes_to_test[1,] #extract the first row to check the function will work with ddpply. We will do it for each row (id row). Using row id we are sure that we are selecting only one row, each combination of genomic interval and transcription factor
test6_function = function(test6_row){

    #extract row index and gene names
    selected_row_id = test6_row$row_index
    selected_known_bat_gene = test6_row$sample_known_bat_genes
    selected_unknown_bat_gene = test6_row$sample_unknown_bat_genes
    selected_random_gene = test6_row$sample_random_genes

    #extract the distance between the known and unknown bat gene
    known_unknown_bat_distance = all_pairs_subset[which(all_pairs_subset$gen1 == as.vector(selected_known_bat_gene) & all_pairs_subset$gen2 == as.vector(selected_unknown_bat_gene)),]$biological_distance

    #extract the distance between the known BAT gene and the random gene
    random_known_bat_distance = all_pairs_subset[which(all_pairs_subset$gen1 == as.vector(selected_known_bat_gene) & all_pairs_subset$gen2 == as.vector(selected_random_gene)),]$biological_distance

    #test if the distance to random is smaller than to unknown BAT gene
    test6 = random_known_bat_distance <= known_unknown_bat_distance

    #bind and return
    return(cbind.data.frame(selected_row_id, selected_known_bat_gene, selected_unknown_bat_gene, selected_random_gene, known_unknown_bat_distance, random_known_bat_distance, test6))
}

#apply the function to raw_tbfs_data_subset, splitting the data.frame for each row_index
test6_final = ddply(.data=genes_to_test, .variables="row_index", .fun=test6_function, .inform=TRUE, .parallel=FALSE, .paropts=NULL)
    #".inform=TRUE" generates and shows the errors. This increases the computation time, BUT is very useful to detect problems in your analyses.
    #".parallel" to paralelize with foreach. 
    #".paropts" is used to indicate additional arguments in for each, specially interesting for using the .export and .packages arguments to supply them so that all cluster nodes have the correct environment set up for computing. 
        #ADD PACKAGES USED INSIDE THE FUNCTION

#THE SCRIPT IS KILLED, TO MUCH MEMORY? MAYBE PARALELZE?

#check the gene ids are correct
if(all(test6_final$row_index == test6_final$selected_row_id)){ #if all TRUE

    #remove the second column with id
    test6_final$selected_row_id <- NULL
} else {

    #if not we have an error
    stop("ERROR!! We have a problem with the gene id in test6_final!!!!")
}

#check we have all the genes
nrow(test6_final) == nrow(genes_to_test) #THIS SHOULD BE THE SAME

#check differences in biological distance
median(test6_final$known_unknown_bat_distance)
median(test6_final$random_known_bat_distance)
sd(test6_final$known_unknown_bat_distance)
sd(test6_final$random_known_bat_distance)

#Calculate p.vale as the proability of a random pair of genes will have lower biological distance than the median of biological distance between all UCP1 connectome genes.
pval = prop.test(x=length(which(test6_final$test6 == TRUE)), n=nrow(test6_final)) #https://stats.stackexchange.com/questions/167164/test-of-equal-proportions-with-zero-successes
pval

#test with a wilcoxon test.
wilcox.test(test6_final$known_unknown_bat_distance, test6_final$random_known_bat_distance, alternative="less") #not sure if we can do this.


#save the results
save.image("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/open_projects/connectome/results/results_2020/rdata/test6_ucp.RData")
require(plyr)
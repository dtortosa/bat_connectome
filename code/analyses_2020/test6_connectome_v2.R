##############################################################################################
####################################### TEST 6 ############################################### 
##############################################################################################

#Test that the Biological distance between known BAT genes and unknown BAT genes inside the connectome are smaller than that of known BAT and random genes.

#we take two BAT genes, we take the median distance between known and non-known genes within the connectome and compare them with the distance between known genes and random sets of genes



###################################################
###### CHANGES RESPECT TO PREVIOUS VERSIONS #######
###################################################

#Respect to version 1: In the previous version I made gene to gene comparisons, but it is too slow, in addition, it does not make sense to make a sampling of the genes within the connectome, because we have a very low number. You should compare each one with random sets, or the median.



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
summary(ucp1_conn) 

#load again the ucp1 connectome again but downloaded in 2020
#Today (19/06/2020), I have downloaded again the UCP1 connectome (the biological distance of human genes to UCP1) and checked that this is exactly similar to the file downloaded 2 years ago and loaded here. 
ucp1_2020 = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/open_projects/connectome/data/human_connectome/UCP1_downloaded_2020.txt", sep="\t", header=TRUE)
identical(ucp1_conn, ucp1_2020) #The files are exactly the same. 



#############################################
######## READ THE BAT INFORMATION ###########
#############################################

#load the information about BAT relationships
bat_relationship = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/open_projects/connectome/results/connectome_results/tables/appendix_S1_ordered.csv", sep=",", header=TRUE)
str(bat_relationship)
summary(bat_relationship)
head(bat_relationship) #we take the supplementary file

#load the new file with the bat relationships obtained in 2020 from the original file of Jose, to check we are using the correct file
bat_relationship_2 = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/open_projects/connectome/results/connectome_results/tables/bat_relationship_check_2020.csv", sep=",", header=TRUE)
#set as zero those cases with 2. I think these are cases not associated but that interact with other proteins that interact with UCP1
bat_relationship_2[which(bat_relationship_2$BAT.relationship == 2),]$BAT.relationship <- 0
#set the rows in alphabetic order according to gene name
bat_relationship_2 = bat_relationship_2[order(bat_relationship_2$Genes),]
#reset the row names
row.names(bat_relationship_2) <- 1:nrow(bat_relationship_2)
#add a new level to the Genes factors
bat_relationship_2$Genes = factor(bat_relationship_2$Genes, c(levels(bat_relationship_2$Genes), "NRIP1"))
#set NRIP1 o RIP140 as NRIP1
bat_relationship_2[which(bat_relationship_2$Genes == "NRIP1 o RIP140"),]$Genes <- "NRIP1"
#remove the not used levels
bat_relationship_2$Genes = droplevels(bat_relationship_2$Genes)

#check
summary(as.vector(bat_relationship$Genes) == as.vector(bat_relationship_2$Genes))
summary(bat_relationship$BAT.relationship == bat_relationship_2$BAT.relationship) #we have the same genes and relationships



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

#set the seed to have reproducible results
set.seed(45788)

#select the distances between known and non-known bat genes.
known_unknown_BAT_distance = all_pairs[which(all_pairs$gen1 %in% known_bat_genes & all_pairs$gen2 %in% unknown_bat_genes),]
#check
summary(known_unknown_BAT_distance$gen1 %in% known_bat_genes)
summary(known_unknown_BAT_distance$gen2 %in% unknown_bat_genes)

##independently to the other, you got the same biological distance
all_pairs[which(all_pairs$gen1 == "UCP1" & all_pairs$gen2 == "FTO"),]
all_pairs[which(all_pairs$gen1 == "FTO" & all_pairs$gen2 == "UCP1"),] #the same biological distance


##first check in more detail
#check that the median distances between known-BAT and candidates genes is the same if you change the order of the genes
median(all_pairs[which(all_pairs$gen2 %in% known_bat_genes & all_pairs$gen1 %in% unknown_bat_genes),]$biological_distance) == median(known_unknown_BAT_distance$biological_distance) #exactly the same distance


##second check in more detail
#select 100 random rows from all pairs
data_for_check_gen1_gen2 = all_pairs[sample(1:nrow(all_pairs), 100),]

#create an empty vector to save results of the check
check_gen1_gen2 = NULL

#fuor each row of data_for_check_gen1_gen2
for(i in 1:nrow(data_for_check_gen1_gen2)){

	#select the [i] row
	selected_row = data_for_check_gen1_gen2[i,]

	#select the row of all pairs that have the gene names inverted, gen1 is gen2 in data_for_check_gen1_gen2 and gen2 is gen1 in data_for_check_gen1_gen2
	opposite_row = all_pairs[which(all_pairs$gen1 == as.vector(selected_row$gen2) & all_pairs$gen2 == as.vector(selected_row$gen1)),] 

	#check and save if both distances are the same
	check_gen1_gen2 = append(check_gen1_gen2, selected_row$biological_distance == opposite_row$biological_distance)
}

#all are the same
check_gen1_gen2 #we assume that a distance from UCP1 to UCP3 is the same than UCP3 to UCP1, therefore, we can subset the all_pairs dataframe with all distances. We will select the all the distances from the known-BAT genes (gen1) to any gene that is not a BAT-candidate. The results should be the same if we took known-BAT genes as gene_2


#subset all pairs to reduce the size
all_pairs_subset = all_pairs[which(all_pairs$gen1 %in% known_bat_genes & !all_pairs$gen2 %in% unknown_bat_genes) ,which(colnames(all_pairs) %in% c("biological_distance"))] 
    #we only want the rows for which the gen1 is a known-bat gene and the gene2 is NOT an unknown bat gene, becuase we want only the random. We test known-BAT with random genes, so we need biological distance columns. We only have now random genes in the second column and known-bat genes in the first one, we can just take random rows.
    #IMPORTANT NOTE. In gen2 we can have known-BAT genes, but this is not problematic. If we have in a random set of distances, distances between known-BAT genes, this would decrease the median, decreasing the power. Even though we not have ANY case where the median random distance was higher than the known-candidate distance (see below). In addition, we are using the median which is less sensitive to outliers, like can be a few known-known distances within a subset of 65000 random distances. If you want to change it you have to add a new conditions "& !all_pairs$gen2 %in% known_bat_genes"
str(all_pairs_subset)

#remove all pairs
rm(all_pairs)

#open an empty vector to save the results
test6 = NULL

#for each iteration
for(i in 1:1000000){

    #randomly select the rows of all_pairs_subset that will be compared with the distance of known - unknown BAT genes
    random_distance_rows = sample(1:length(all_pairs_subset), nrow(known_unknown_BAT_distance), replace=FALSE) 
    	#we need the same size than known_unknown_BAT_distance, for having two comparable datasets
    	#we need replace=FALSE, because each iteration, i want to have different 6527 distances, not repeated cases. You can have repeated values between iterations, but this is not a problem, because we are selecting 6000 from 1 million, so the odds are not very high. 

    #select the rows from all_pairs_subset
    random_known_BAT_distance = all_pairs_subset[random_distance_rows]

    #check if the median distance between random and known-BAT is smaller or equal than the median distance of known - unknown BAT genes
    test6 = append(test6, median(random_known_BAT_distance) <= median(known_unknown_BAT_distance$biological_distance))
}


#Calculate p.vale as the proability of a random pair of genes will have lower biological distance than the median of biological distance between all UCP1 connectome genes.
pval = prop.test(x=length(which(test6 == TRUE)), n=length(test6)) #https://stats.stackexchange.com/questions/167164/test-of-equal-proportions-with-zero-successes
pval

#save the results
save.image("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/open_projects/connectome/results/results_2020/rdata/test6_ucp_v2.RData")
require(plyr)
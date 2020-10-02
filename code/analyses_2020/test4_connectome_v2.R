####################################
############# TEST 4 ###############
####################################

#Scripts for testing if number of the known BAT genes include in BAT connectome is higher than expected by chance. 


#HEMOS QUITADO EL TEST 4. 
	#esto sirve para testar que por azar no se explica haber conseguido justo esos 60 genes, que no es lo ideal, pero no tenemos una lista completa de bat genes a lo largo del genoma.
		#por eso hemos hecho el test 7, testando que los known BAT genes estan mas cerca entre respecto de random genes. Esto lo sugirió Yuval.
	#Sobre los problemas del Test 4: Si tener esos 60 known-BAT genes inside the BAT connectome is expected by chance.
		#We are checking in that way if we can get these 60 genes by chance just selecting a random set of 167 genes. This is not the best test. The best test would be check if random sets of genes have 60 or more known BAT genes (included or not in the connectome), but this would entail to curate all the human genes. We have the good and strong test comparing the distance between known-BAT and candidate genes with known-BAT and random genes. This is test 6 and compare the distance within the connectome but specifying between the two relevant groups, known and candidates. This is an upgrade of the initial test comparing distance inside and outside of the connectome. 



###################################################
###### CHANGES RESPECT TO PREVIOUS VERSIONS #######
###################################################

#Respect to version 1: In that version we calculate the proportion of random sets of 60 genes that have 60 known BAT genes. This is not correct. As Yuval pointed out, we need to select random sets of 167 genes (size of the BAT connectome) and check how many had our 60 known BAT genes inside. In that way we are checking if we can get these 60 genes by chance just selecting a random set of 167 genes. This is not the best test. The best test would be check if random sets of genes have 60 or more known BAT genes (included or not in the connectome), but this would entail to curate all the human genes. We have the good and strong test comparing the distance between known-BAT and candidate genes with known-BAT and random genes. This is test 6 and compare the distance within the connectome but specifying between the two relevant groups, known and candidates. This is an upgrade of the initial test comparing distance inside and outside of the connectome. 
 


#############################################
######## READ THE UCP1 CONNECTOME ###########
#############################################
#load the connectome
ucp1_conn = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/open_projects/connectome/data/human_connectome/UCP1.txt", sep="\t", header=T)
str(ucp1_conn)
head(ucp1_conn)

#create a vector with all codifican human genes ordered alphabetically
all_genes = sort(ucp1_conn$Target)

#extract top 1% of UCP1 connectome 
BAT_connectome = ucp1_conn[which(ucp1_conn$Target_in_source_P.value.percentile<0.01),]
length(BAT_connectome$Target)


### read data of function related with BAT of BAT connectome
##read the old file
func_1 = read.csv("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/open_projects/connectome/code/connectome/tree/relation_with_BAT.csv", header=T)
#set as zero those cases with 2. I think these are cases not associated but that interact with other proteins that interact with UCP1
func_1[which(func_1$FUN == 2),]$FUN <- 0
#convert FUN to factor
func_1$FUN = as.factor(func_1$FUN)
#take a look
str(func_1)

##read the new file used in the test 6 in 2020
func_2 = read.csv("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/open_projects/connectome/results/connectome_results/tables/appendix_S1_ordered.csv", header=T)
#change colnames
colnames(func_2)[which(colnames(func_2) == "Genes")] <- "GENES" 
colnames(func_2)[which(colnames(func_2) == "BAT.relationship")] <- "FUN"
#convert FUN to factor
func_2$FUN = as.factor(func_2$FUN) 
#take a look
str(func_2)

##merge both files and check they have the same FUN data
#merge
merge_test = merge(func_1, func_2, by="GENES")
#check FUN data
summary(merge_test$FUN.x == merge_test$FUN.y)
    #They are the same, so we will continue using the same file used in the previous script version

##now check that the genes included in the BAT connectome according to these file are the same than those included in the connectome obtained directly from UCP1.txt
summary(as.vector(merge_test$GENES) %in% as.vector(BAT_connectome$Target))
summary(as.vector(BAT_connectome$Target) %in% as.vector(merge_test$GENES))

#check the same for the list from all_genes
summary(BAT_connectome$Target %in% all_genes[which(all_genes %in% func_1$GENES)])
summary(all_genes[which(all_genes %in% func_1$GENES)] %in% BAT_connectome$Target)


### drop UCP1 because UCP1 was used as core genes, thus of course will be related with BAT, we used because of this! Therefore, we will have a list of 167 genes related with according to the bibliography
BAT_connectome = BAT_connectome[-which(BAT_connectome$Target == "UCP1"),]
#take a lool
str(BAT_connectome)
head(BAT_connectome)

### select those genes of the BAT connectome that are associated with BAT according evidence of 2018
BAT_known = BAT_connectome[which(BAT_connectome$Target %in% func_1[which(func_1$FUN == 1),]$GENES),]
#check
summary(BAT_known$Target %in% func_1[which(func_1$FUN == 1),]$GENES)


### loop to run the simulations
#set the seed
seed = 345767
set.seed(seed)

#set the number of iterations
n_iter = 1000000

#run the loop
results = NULL
for(i in 1:n_iter){

    ### random selection of as many genes as those included in the BAT connectome
    random_genes = all_genes[sample(1:length(all_genes), nrow(BAT_connectome))]
    #random_genes

    ### from random genes select those that are known to be related to the BAT and are included in the BAT connectome 
    random_known_bat = random_genes[which(random_genes %in% BAT_known$Target)]

    ### check if ALL the genes that are known-BAT genes and are included in the connectome are also included in the random connectome. We are checking in that way if we can get these 60 genes by chance just selecting a random set of 167 genes. This is not the best test. The best test would be check if random sets of genes have 60 or more known BAT genes (included or not in the connectome), but this would entail to curate all the human genes. We have the good and strong test comparing the distance between known-BAT and candidate genes with known-BAT and random genes. This is test 6 and compare the distance within the connectome but specifying between the two relevant groups, known and candidates. This is an upgrade of the initial test comparing distance inside and outside of the connectome. 
    results = append(results, length(random_known_bat) == nrow(BAT_known))
}

#Calculate p.vale as the probability of a random set of genes have equal or higher number of BAT related genes
pval = prop.test(x=length(which(results == TRUE)), n=length(results)) #https://stats.stackexchange.com/questions/167164/test-of-equal-proportions-with-zero-successes

#print pval 
pval

#save results
write.table(results, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/open_projects/connectome/results/results_2020/tests/test_4/results_test_4_v2.txt", col.names = FALSE, row.names = FALSE, sep="\t")

#save the results
save.image("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/open_projects/connectome/results/results_2020/rdata/test_4_v2.RData")
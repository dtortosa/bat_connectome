#Test 2: Biological distance of random set of of genes is lower than the median biological distance than the UCP1 connectome

####load data of all the pairs across genome
load("/home/dsalazar/ucp_connectome/data/connectome_data_loaded.RData")

#############################################
######## READ THE UCP1 CONNECTOME ###########
#############################################
ucp1_conn = read.table("/home/dsalazar/ucp_connectome/data/UCP1.txt", sep="\t", header=T)
str(ucp1_conn)
head(ucp1_conn)


#############################################
######## SELECT TOP 1% BY P.VALUE ###########
#############################################

##Select the top 1% of all human genes by p-value, only genes close to the core by a p.value lower than 0.01. You have to use the core p.value, which is the p.value about the distance between the core of your connectome and the target gen. The next column (target to core), is the p.value from the target to the core gene but in the connectome of the target gene. For example, Target_in_source_P.value.percentile. of UCP2 in the connectome of UCP1 is 0.00006, whilst the Source_in_target_P.value.percentile. is 0.00012, exactly the Target_in_source_P.value.percentile. of UCP1 in the connectome of UCP2. 
subset = ucp1_conn[ucp1_conn$Target_in_source_P.value.percentile<0.01,]
subset$Target
length(subset$Target)

##########################
######## TEST 2 ##########
##########################

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

####select those rows form all_pairs in which the pair is included in the UCP1 connectome
pairs_ucp1_connec = all_pairs[which(all_pairs$pairs %in% ucp_conn_pairs),]

####calculate median of biological distance between UCP1 connectome genes
ucp1_median_dist = median(pairs_ucp1_connec$biological_distance)

###create a vector with all codifican human genes ordered alphabetically
all_genes = sort(ucp1_conn$Target)

###loop for testing if the proportion of simulated random genes set with a median biological distance smaller than for the UCP1 connectome genes


#function of the loop
test2_loop = function(n_iter, seed){

    #create a file for resutls
    file.create(paste("/home/dsalazar/ucp_connectome/results/test_2/results_", seed, ".txt", sep=""))

    #set the seed
    set.seed(seed)

    #loop
    for (i in 1:n_iter){
    
        ###random selection of as many genes as those included in UCP1 connectome
        random_genes = all_genes[sample(1:length(all_genes), nrow(subset))]
    
        #####Loop for createing a vector with all possible pairs of genes included en in the top 1% of UCP1 connectome.
        random_genes_pairs = NULL
        for(i in 1:length(random_genes)){ #for each row of subset
    
            #select the row [i]
            selected_gen = random_genes[i] 
    
            #if the it is not the final row
            if(!i==length(random_genes)){ 
    
                #select all the rows following the [i] row, in this way each pair is only written one time 
                subset_2 = random_genes[(i+1):168]  
    
                #for each gen of subset2    
                for(j in 1:length(subset_2)){ 
    
                    #select the gen [j]
                    selected_gen_2 = subset_2[j] 
    
                    #create a string with the [i] gen and the [j] gen.
                    random_genes_pairs = append(random_genes_pairs, paste(selected_gen, selected_gen_2, sep="_")) 
                }
            }
        }
    
        ####select those rows form all_pairs in which the pair is included in the random genes
        pairs_random_connec = all_pairs[which(all_pairs$pairs %in% random_genes_pairs),]
    
        ####calculate median of biological distance between random genes genes
        random_median_dist = median(pairs_random_connec$biological_distance)
    
        #### test of distance between genes of the random random genes are lower than the the median distance of genes of the UCP1 connectome 
        if(random_median_dist < ucp1_median_dist){
            result_test2 = "YES"
        } else {
            result_test2 = "NO"
        } 

        #write the result of [i] iteration
        write.table(result_test2, paste("/home/dsalazar/ucp_connectome/results/test_2/results_", seed, ".txt", sep=""), sep = "\n", append=TRUE, col.names = FALSE, row.names = FALSE) #"\n" for obtain each results in the same column.
    } 
}

#Paralelize the process
require(foreach) #for repeat a process several times 
require(doParallel) #for parallel

# set up cluster
clust <- makeCluster(3, type="FORK") 
registerDoParallel(clust)

#create a vector with seed of each process
seed = 1:3

##run each process separately
foreach(seed = seed) %dopar% { 
    test2_loop(seed = seed, n_iter = 3334)
} 

#stop the cluster 
stopCluster(clust)

#eliminate the bid file with all pairs of genes
rm(all_pairs)
 
#save the results
save.image("/home/dsalazar/ucp_connectome/results/rdata/test2_ucp.RData")


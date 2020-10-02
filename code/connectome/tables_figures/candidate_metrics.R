#La idea de este codigo era ver cual es el gen más cercano para cada uno de los candidatos cogiendo la distanca biológica más corta. El problema es que no cuadra con la filogenia, y dado que en esa filogenia se tienen en cuenta las distnacias del gen candaidato, pero tambén todas las distancias del que "debería" estar más cerca (con el candidato y muchos otros), creo que es más fina la filogenia. A parte, yo la hice indicando una lista de genes (los 168 del conectoma), no toqué nada más, el resto lo hizo el codigo de Yuval. 

####load data of all the pairs across genome
load("/Users/diegosalazar/phd_big_documents/UCP_connectome/connectome_data_loaded.RData")

#############################################
######## READ THE UCP1 CONNECTOME ###########
#############################################
ucp1_conn = read.table("/Users/diegosalazar/Google Drive/science/open_projects/connectome/data/human_connectome/UCP1.txt", sep="\t", header=T)
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

####select those rows form all_pairs in which the pair is included in the UCP1 connectome
pairs_ucp1_connec = all_pairs[which(all_pairs$pairs %in% ucp_conn_pairs),]


### read data of function related with BAT of all UCP1 connectome genes
func = read.csv("/Users/diegosalazar/Google Drive/science/open_projects/connectome/code/connectome/tree/relation_with_BAT.csv", header=T)
str(func)
func$FUN = as.factor(func$FUN)

candidate = "LPAR2"

distances = pairs_ucp1_connec[which(pairs_ucp1_connec$gen1 == candidate & pairs_ucp1_connec$gen2 %in% func[which(func$FUN == 1),]$GENE | pairs_ucp1_connec$gen2 == candidate & pairs_ucp1_connec$gen1 %in% func[which(func$FUN == 1),]$GENE),]

pval = NULL
for(i in 1:nrow(distances)){
    selected_row = distances[i,]
    if(selected_row$gen1 == candidate){
        pval = append(pval, selected_row$Target_in_source_P.value.percentile)
    } else {
        pval = append(pval, selected_row$Source_in_target_P.value.percentile)
    }
}

final_distances = cbind.data.frame(distances, pval)

final_distances[which(final_distances$pval == min(final_distances$pval)),]



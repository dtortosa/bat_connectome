#############################################
######## READ ALL PAIRS DATA ################
#############################################
all_pairs = read.table("/Users/diegosalazar/phd_big_documents/UCP_connectome/combined_all_dist_rank_p_10.00.txt", sep="\t", header=T)
colnames(all_pairs) <- c("gen1", "gen2", "biological_distance", "connectome_rank", "Target_in_source_P.value.percentile", "Source_in_target_P.value.percentile")
str(all_pairs)
head(all_pairs)

####save data of all the pairs across genome
save.image("/Users/diegosalazar/phd_big_documents/UCP_connectome/connectome_data_loaded.RData")

####create a colum with pairs between genes
all_pairs$pairs = paste(all_pairs$gen1, all_pairs$gen2, sep="_")

####select only pairs and biological distance variables and overwrite all_pairs (less size)
all_pairs = all_pairs[, which(colnames(all_pairs) %in% c("pairs", "biological_distance"))]

####save data of all the pairs across genome
save.image("/Users/diegosalazar/phd_big_documents/UCP_connectome/connectome_data_loaded_with_pairs.RData")

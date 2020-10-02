#############################################
######## READ ALL PAIRS DATA ################
#############################################
load("/home/dsalazar/ucp_connectome/data/connectome_data_loaded.RData")

####create a colum with pairs between genes
all_pairs$pairs = paste(all_pairs$gen1, all_pairs$gen2, sep="_")

####select only pairs and biological distance variables
all_pairs = all_pairs[, which(colnames(all_pairs) %in% c("pairs", "biological_distance"))]

####save data of all the pairs across genome
save.image("/home/dsalazar/ucp_connectome/data/connectome_data_loaded.RData")

#code for reducing thr size of connectome_data_loaded.RData dropping all columns except pairs and biological_distance

####load data of all the pairs across genome by menad of a RData (faster)
load("/SCRATCH/UGR002/dsalazar/ucp1_connectome/data/connectome_data_loaded.RData")

####select only pairs and biological distance variables and overwrite all_pairs (less size)
all_pairs = all_pairs[, which(colnames(all_pairs) %in% c("pairs", "biological_distance"))]

####save data of all the pairs across genome
save.image("/SCRATCH/UGR002/dsalazar/ucp1_connectome/data/connectome_data_loaded_with_pairs.RData")

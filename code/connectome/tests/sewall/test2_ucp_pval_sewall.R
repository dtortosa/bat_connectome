#code for extract p.values 

#select path of results of all parallel process
result_paths = list.files("/home/dsalazar/ucp_connectome/results/test_2", pattern="results_", full.names = TRUE)

#read each result
result_test2 = as.data.frame(NA)
colnames(result_test2)[1] <- "V1"
for(i in 1:length(result_paths)){
    path_selected = result_paths[i]
    result = read.table(path_selected, header=F, sep="\t")
    result_test2 = rbind(result_test2, result)
}
result_test2 = result_test2[-1,]

#Calculate p.vale as the proability of a random pair of genes will have lower biological distance than the median of biological distance between all UCP1 connectome genes.
pval = prop.test(x=length(which(result_test2 == "YES")), n=length(result_test2)) #https://stats.stackexchange.com/questions/167164/test-of-equal-proportions-with-zero-successes
#eliminate the bid file with all pairs of genes

rm(all_pairs)

#save the results
save.image("/home/dsalazar/ucp_connectome/results/test2_ucp_pval.RData")
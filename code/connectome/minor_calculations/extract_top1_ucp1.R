########Code for working in the BAT connectome#########

####We only will use the connectome of UCP1. Some proteins are under discussion yet (e.g. UCP2, 3), whilst other proteins, such us PPARGC1A, are implicated in heat production in BAT, but also in heat production outside of this tissue. The advice of Yuval is "I would usually start with the genes that have the best experimental evidence so far in being involved in the biological process that you are interested in. In that case UCP1 could be a good starting point. If you want to test computationally whether the other genes are good candidates, you can test their p-value (or biological distance) with UCP1, and in that case I see that all have a very low p-value. However, like you mentioned it may be that the tissue of expression of the other genes may not be relevant, and in that case the prediction would be less reliable. So I would start with the connectome UCP1 only (or if there are other genes with strong experimental evidence), capturing (say) all genes with p<0.01, and then you can further filter by excluding the genes that are not expressed in the relevant tissues (you can use resources such as bioGPS or GTEx for that)." 
#Jonantan and Jose support the idea of using only UCP1 as core. 
#You have to make a filter by p.value and then by expression in the tissue (only genes expressed in the tissue)

#############################################
######## READ THE UCP1 CONNECTOME ###########
#############################################
ucp1_conn = read.table("/Users/diegosalazar/Google Drive/science/open_projects/bat_geno_expression/data/human_connectome/UCP1.txt", sep="\t", header=T)
str(ucp1_conn)
head(ucp1_conn)


#############################################
######## SELECT TOP 1% BY P.VALUE ###########
#############################################

##Select the top 1% of all human genes by p-value, only genes close to the core by a p.value lower than 0.01. You have to use the core p.value, which is the p.value about the distance between the core of your connectome and the target gen. The next column (target to core), is the p.value from the target to the core gene but in the connectome of the target gene. For example, Target_in_source_P.value.percentile. of UCP2 in the connectome of UCP1 is 0.00006, whilst the Source_in_target_P.value.percentile. is 0.00012, exactly the Target_in_source_P.value.percentile. of UCP1 in the connectome of UCP2. 
subset = ucp1_conn[ucp1_conn$Target_in_source_P.value.percentile<0.01,]
subset$Target
length(subset$Target)

#save the list 
write.table(subset$Target, "/Users/diegosalazar/Google Drive/science/open_projects/bat_geno_expression/results/connectome_results/UCP1_connectome_less_0.01.txt", row.names=FALSE, col.names = FALSE)

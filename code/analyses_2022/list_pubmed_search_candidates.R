library(easyPubMed)
	library(R.utils)


#############################################
######## READ THE BAT INFORMATION ###########
#############################################

#load the information about BAT relationships
bat_relationship = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/connectome_results/tables/appendix_S1_ordered.csv", sep=",", header=TRUE)
str(bat_relationship)
summary(bat_relationship)
head(bat_relationship) #we take the supplementary file

#load the new file with the bat relationships obtained in 2020 from the original file of Jose, to check we are using the correct file
bat_relationship_2 = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/connectome_results/tables/bat_relationship_check_2020.csv", sep=",", header=TRUE)
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

#select those genes that are NOT known to be associated with BAT within the connectome
unknown_bat_genes = bat_relationship[which(bat_relationship$BAT.relationship == 0),]$Genes



list_abstracts = list()
for(i in 1:length(unknown_bat_genes)){

	selected_unknown_bat_gene = unknown_bat_genes[i]

	my_query <- paste('(((("Adipose Tissue, Brown"[Mesh] OR "Brown Fat" OR "Brown adipose tissue"))) OR (("Adipose tissue, beige"[Mesh] OR "beige adipose tissue" OR "Brite fat" OR "beige fat")))AND(', selected_unknown_bat_gene, ')', sep="")
	my_entrez_id <- get_pubmed_ids(my_query)
	#my_abstracts_txt <- fetch_pubmed_data(my_entrez_id, format = "abstract")
	#my_abstracts_txt
	
	#https://cran.r-project.org/web/packages/easyPubMed/vignettes/getting_started_with_easyPubMed.html

	my_abstracts_xml = withTimeout(fetch_pubmed_data(pubmed_id_list = my_entrez_id), timeout=2, onTimeout="warning")
		#https://stackoverflow.com/questions/31462416/r-set-execution-time-limit-in-loop

	my_titles <- custom_grep(my_abstracts_xml, "ArticleTitle", "char")
	
	# use gsub to remove the tag, also trim long titles
	#TTM <- nchar(my_titles) > 75
	#my_titles[TTM] <- paste(substr(my_titles[TTM], 1, 70), "...", sep = "")
	
	# Print as a data.frame (use kable)
	head(my_titles)

	if(!grepl("reached", my_abstracts_xml)){
		list_abstracts[[i]] <- my_titles
	} else {
		list_abstracts[[i]] <- NA	
	}
}

names(list_abstracts) = unknown_bat_genes




save.image("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/results_2022/list_pubmed_search_candidates.RData")


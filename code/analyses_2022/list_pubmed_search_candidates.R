#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



######################################################################################
####################### PUBMED SEARCH CANDIDATES #####################################
######################################################################################

#extract the title of papers obtained for a search about each of the 107 BAT candidates and BAT functioning



#################################################################
####################### REQUIRED PACKAGES #######################
#################################################################

library(easyPubMed) #to connect to pubmed: 	
	#https://cran.r-project.org/web/packages/easyPubMed/vignettes/getting_started_with_easyPubMed.html
library(R.utils) #to stop processes taking too much time



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



#####################################
######## MAKE THE SEARCH  ###########
#####################################

#open an empty list for the abstract
list_abstracts = list()

#for each unknown BAT gene, i.e., candidate
for(i in 1:length(unknown_bat_genes)){

	#select the [i] candidate
	selected_unknown_bat_gene = unknown_bat_genes[i]

	#use the exact same query used my Jose when looking for gen-BAT association
	my_query <- paste('(((("Adipose Tissue, Brown"[Mesh] OR "Brown Fat" OR "Brown adipose tissue"))) OR (("Adipose tissue, beige"[Mesh] OR "beige adipose tissue" OR "Brite fat" OR "beige fat")))AND(', selected_unknown_bat_gene, ')', sep="")
	
	#make the query
	my_entrez_id <- get_pubmed_ids(my_query)
		#Query PubMed (Entrez) in a simple way via the PubMed API eSearch function. Calling this function results in posting the query results on the PubMed History Server. This allows later access to the resulting data via the fetch_pubmed_data() function, or other easyPubMed functions.
	
	#extract the titles of the papers with fetch_pubmed_data
	my_titles_xml = withTimeout(fetch_pubmed_data(pubmed_id_list = my_entrez_id), timeout=5, onTimeout="warning")
		#we make an envelop with withTimeout in order to avoid too much time looking. I have detected that the function gets stuck if there are no results. We limit the execution time to 2 seconds. If the limit is surpass, a warning is generated and save in the corresponding object. 
			#https://stackoverflow.com/questions/31462416/r-set-execution-time-limit-in-loop
		#you can also use fetch_pubmed_data to load the abstracts
			#my_abstracts_txt <- fetch_pubmed_data(my_entrez_id, format = "abstract")

	#extract the titles from the xml file
	my_titles <- custom_grep(my_titles_xml, tag="ArticleTitle", format="char")
		#Extract text form a string containing XML or HTML tags. Text included between tags of interest will be returned. If multiple tagged substrings are found, they will be returned as different elements of a list or character vector.

	#use gsub to remove the tag, also trim long titles
	#TTM <- nchar(my_titles) > 75
	#my_titles[TTM] <- paste(substr(my_titles[TTM], 1, 70), "...", sep = "")
	
	#if we do not have the warning saved in my_titles_xml
	if(!grepl("reached", my_titles_xml)){

		#save the titles as the [[i]] element of the list
		list_abstracts[[i]] <- my_titles
	} else { #if not, and then, not title was obtained

		#save NA in the [[i]] position
		list_abstracts[[i]] <- NA	
	}

	#set the gene name
	names(list_abstracts)[[i]] = selected_unknown_bat_gene
}

#check
length(list_abstracts) == length(unknown_bat_genes)
!FALSE %in% c(names(list_abstracts) == unknown_bat_genes)


save.image("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/results_2022/list_pubmed_search_candidates_2.RData")
#load("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/results_2022/list_pubmed_search_candidates.RData")
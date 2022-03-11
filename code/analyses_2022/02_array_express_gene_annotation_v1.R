#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



####################################################################################
########################## PREPROCESING GENE EXPRESSSION ###########################
####################################################################################

#Script for obtaining gene expression from ArrayExpress processed data

#Part of the code of this script comes from
	#a tutorial of ArrayExpress
		#http://www.bioconductor.org/packages/release/bioc/vignettes/ArrayExpress/inst/doc/ArrayExpress.pdf
	
	#a book chapter about using CEL files to get genes different between populations, tissues..
		#From CEL Files to Annotated Lists of Interesting Genes



#################################################################
####################### REQUIRED PACKAGES #######################
#################################################################

require(ArrayExpress) #for loading data from ArrayExpress
require(arrayQualityMetrics) #for normalizing expression data in each dataset
require(dplyr) #for doing some data operations
require(plyr) #for llply
require(foreach) #for parallel
require(doParallel) #for parallel



########################################################
####################### STARTING #######################
########################################################

#set working directory to results, because ArrayExpress save stuff during the analyses automatically
wd="/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/results_2022"
setwd(wd)

#completely remove the previous folder with raw data
system(paste("rm -rf array_express_gene_expression", sep=""))

#create a directory to save the results and move to it
system(paste("mkdir -p array_express_gene_expression", sep=""))
	#use mkdir -p because if you don't and run again, you will get an error that it exists. "p" flags is for "no error if existing, make parent directories as needed"

#load the RDS file with the information of all BAT datasets included in the analysis
bat_datasets_info = readRDS("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/results_2022/array_express_raw/info_datasets.Rds")

#extract the name of the datasets included
ids_bat_studies = as.vector(unlist(sapply(X=bat_datasets_info, "[", "dataset")))



################################################################################
####################### WRITE FUNCTION FOR FULL ANALYSIS #######################
################################################################################


##function for doing all the analyses across the datasets of interest
#bat_studies_considered=ids_bat_studies[which(ids_bat_studies %in% c("E-GEOD-27657", "E-GEOD-54280"))] #for debugging
full_analysis = function(bat_studies_considered){

	#print start
	print("###############################################")
	print(paste("STARTING: ", paste(bat_studies_considered, collapse=" | "), sep=""))
	print("###############################################")



	##################################################################
	####################### EXPRESSION BY GENE #######################
	##################################################################

	##function for calculating expression per gene
	#selected_ids_bat_studies=bat_studies_considered[1] #for debugging
	extract_expression_data = function(selected_ids_bat_studies){
	
		##starting
		#make a directory for the selected dataset
		system(paste("mkdir -p array_express_gene_expression/", selected_ids_bat_studies, sep=""))
			#use mkdir -p because if you don't and run again, you will get an error that it exists. "p" flags is for "no error if existing, make parent directories as needed"
	
		#see one specific dataset
		print("###############################################")
		print(paste("STARTING ", selected_ids_bat_studies, ": ", sep=""))
		print("###############################################")
	
		#load the AE set filtered
		AEsetnorm_filter = readRDS(paste("array_express_raw/", selected_ids_bat_studies, "/", selected_ids_bat_studies, "_expression_filter_1.Rds", sep=""))
	
	
		##prepare gene expression data
		#get a data.frame with rows representing probes and columns representing arrays. To know which columns go with that tissue or experimental treatment, we can rely on the phenoData information inherited from AEset_human.
		expression_matrix = data.frame(exprs(AEsetnorm_filter))
		#head(expression_matrix)
			#access the expression and error measurements of assay data stored in an object derived from the ‘eSet-class’.
		#head(pData(AEsetnorm_filter))
	
		#we have a value of expression per sample, so we have to find a way to summarize and get a value across all samples per probe We have filtered the samples to be considered in a previous script, so we can just use the median to summarize the expression across all the remaining samples. These should be only BAT-related samples of humans thanks to the filtering. In the same experiments, we can have BAT sample of males and females of different ages, or brown adipocytes generated in different ways. But thanks to the filtering, we know that the remaining samples belong to experimental groups that show BAT-like features.
		#We will use median like in the rest of the paper in order to summarize.
		
		#calculate the median expression across all samples per prob
		median_expression = as.data.frame(apply(X=expression_matrix, MARGIN=1, FUN=median))
		#set the column name as average expression
		colnames(median_expression) = "average_gene_expression"
		#str(median_expression)
		#check
		print("###############################################")
		print(paste("MERGIN", selected_ids_bat_studies, " OK?: ", sep="")); print(identical(row.names(median_expression), row.names(expression_matrix)))
		print("###############################################")
	
		#merge the result with the original data.frame using the row names 
		expression_matrix = merge(expression_matrix, median_expression, by="row.names")
		#str(expression_matrix)
	
		#select only the row names and the average gene expression
		expression_matrix_average = expression_matrix[,c("Row.names", "average_gene_expression")]
		#change the colname of the probeIDS
		colnames(expression_matrix_average)[which(colnames(expression_matrix_average) == "Row.names")] = "PROBEID"
		#str(expression_matrix_average)
		#head(expression_matrix_average)
	
	
		##gene annotation
		#Traditionally, Affymetrix arrays (the so-called 3’ IVT arrays) were probeset based: a certain fixed group of probes were part of a probeset which represented a certain gene or transcript (note however, that a gene can be represented by multiple probesets). The more recent “Gene” and “Exon” Affymetrix arrays are exon based and hence there are two levels of summarization to get to the gene level. The “probeset” summarization leads to the exon level. The gene / transcript level is given by “transcript clusters”. Hence, the appropriate annotation package for our chip type should be something like hugene10sttranscriptcluster.db.
			#In general, I have avoided annotation packages that have the name probset, like hugene10stprobeset.db.
			#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#74_Old_and_new_%E2%80%9Cprobesets%E2%80%9D_of_Affymetrix_microarrays
	
		#We use the function mapIds from AnnotationDbi to query the gene symbols for the transcript clusters. For each cluster, we added the gene symbol (SYMBOL)
			#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#11_Annotation_of_the_transcript_clusters	
	
		#select the AnnotationDb object
		#we have to convert the probs ID to the gene symbols used by Yuval. In all cases, we select the annotation package that it is not probset. I understand that we do not need probset summarization, which leads to exons, but gene summarization
			#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#74_Old_and_new_%E2%80%9Cprobesets%E2%80%9D_of_Affymetrix_microarrays
		if(AEsetnorm_filter@annotation == "pd.hugene.1.0.st.v1"){
	
			#we select the hugene 10 annotation package
			require(hugene10sttranscriptcluster.db)
			selected_annotation = hugene10sttranscriptcluster.db
				#we have avoided hugene10stprobeset.db
					#http://bioconductor.org/packages/release/BiocViews.html#___AnnotationData
						#search for "hugene"
					#http://bioconductor.org/packages/release/data/annotation/manuals/hugene10sttranscriptcluster.db/man/hugene10sttranscriptcluster.db.pdf
		}
		if(AEsetnorm_filter@annotation == "pd.hg.u133.plus.2"){
			
			#we select the hg 133 annotation package
			require(hgu133plus2.db)
			selected_annotation = hgu133plus2.db
				#we have avoided hgu133plus2probe (see above)
					#http://bioconductor.org/packages/release/BiocViews.html#___AnnotationData
						#search for "133"
					#http://bioconductor.org/packages/release/data/annotation/manuals/hgu133plus2.db/man/hgu133plus2.db.pdf
		}
	
		#using the selected annotation package, extract the gene symbols for the probs we have (rows in the expression data)
		gene_symbols_probs = AnnotationDbi::select(x=selected_annotation, keys=expression_matrix_average$PROBEID, column=c("SYMBOL"), keytype="PROBEID")
			#x: AnnotationDb object
			#keys: The keys to select records from the database.
				#we are using the probes IDs of the selected dataset
			#columns: the columns or kinds of things that can be retrieved from the database.
			#keytype: the keytype that matches the keys used. 
				#We are using the probe ids.
	
		
		##remove NA
		#probes with and without SYMBOL
		probs_without_symbol = which(is.na(gene_symbols_probs$SYMBOL))
		probs_with_symbol = which(!is.na(gene_symbols_probs$SYMBOL))
		
		#see
		print("###############################################")
		print(paste("FROM ", nrow(gene_symbols_probs), " PROBS, ", length(probs_without_symbol), " HAVE NO GENE SYMBOL", sep=""))
		print("###############################################")
	
		#select those probes with gene symbol
		gene_symbols_probs = gene_symbols_probs[probs_with_symbol,]
	
		#check
		print("###############################################")
		print(paste("ALL PROBS WITHOUT GENE SYMBOL WERE REMOVED?:", sep="")); print(length(which(is.na(gene_symbols_probs$SYMBOL))) == 0)
		print("###############################################")
	
	
		##remove probs associated with several genes
			#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#111_Removing_multiple_mappings
		#first put together rows with the same PROBEID as a group
		probes_grouped <- group_by(gene_symbols_probs, PROBEID)
			#Most data operations are done on groups defined by variables. group_by() takes an existing tbl and converts it into a grouped tbl where operations are performed "by group". ungroup() removes grouping.
	
		#then summarize by counting the number of distinct symbols per probe ID
		probes_summarized <- dplyr::summarize(.data=probes_grouped, no_of_matches=n_distinct(SYMBOL))
			#‘summarise()’ creates a new data frame. It will have one (or more) rows for each combination of grouping variables
			#n_distinct:  This is a faster and more concise equivalent of ‘length(unique(x))
			#For each group (i.e., each probe), calculate the number of unique gene symbols associated. This is the number of genes that each probe is associated with. Save the result in no_of_matches.
		#head(probes_summarized)
	
		#select those probes associated with more than 1 gene symbol
		probe_stats <- filter(probes_summarized, no_of_matches > 1)
			#The ‘filter()’ function is used to subset a data frame, retaining all rows that satisfy your conditions. 
	
		#see
		print("###############################################")
		print(paste("HOW MANY PROBES ARE ASSOCIATED WITH SEVERAL GENE SYMBOLS?: ", sep="")); print(nrow(probe_stats))
		print("###############################################")
	
		#remove these probes
		gene_symbols_probs = gene_symbols_probs[which(!gene_symbols_probs$PROBEID %in% probe_stats$PROBEID),]
			#We have probes that map to multiple gene symbols. It is difficult to decide which mapping is “correct”. Therefore, we exclude these probes.
	
		#check
		print("###############################################")
		print(paste("WE HAVE REMOVED ALL PROBES ASSOCIATED WITH MULTIPLE GENES: ", sep="")); print(unique(dplyr::summarize(.data=group_by(gene_symbols_probs, PROBEID), no_of_matches=n_distinct(SYMBOL))$no_of_matches) == 1) #see previous lines for an explanation of these steps
		print("###############################################")
	
	
		##merge annotation and expression data
		#merge by PROBEID
		merged_data = merge(expression_matrix_average, gene_symbols_probs, by="PROBEID") #CHECK THIS
		#str(merged_data)
		#head(merged_data)
		
		#see
		print("###############################################")
		print(paste("SEE SUMMARY: ", sep="")); print(summary(merged_data)) #see previous lines for an explanation of these steps
		print("###############################################")
	
		
		##summarize expression by gene
		#The same gene can be detected by several probes, because the same gene can produce different transcripts. We consider the median of all transcript per gene
	
		#there is some controversy
			#If they are different isoforms, doing an average of two might not be appropriate cause one might no be expressed and that brings gene expression down..
			#Yes, of course, in which case you will have to write some code to check for these situations in which the expression is so low for one probe such that it is negligible.
				#https://www.biostars.org/p/271379/#271588
			#otros dicen the usar el valor maximo
				#https://support.bioconductor.org/p/70133/
			#we are just going to use the usual approach with median, which is less sensitive to outliers.
	
		#see cases
		print("###############################################")
		print(paste("CASES IN WHICH A GENE IS DETECTED BY SEVERAL PROBES: ", sep="")); print(length(which(duplicated(merged_data$SYMBOL))))
		print("###############################################")
	
		#calculate the median expression in each gene
		final_data = aggregate(average_gene_expression ~ SYMBOL, merged_data, median)
			#some people recommend to summarize probes by gene using mean or median. We have used median to summarize across the whole manuscript, so we will stick to it.
				#https://www.biostars.org/p/271379/#271588
				#https://www.biostars.org/p/336130/
	
		#see duplicated cases
		print("###############################################")
		print(paste("CASES IN WHICH A GENE IS DETECTED BY SEVERAL PROBES AFTER SUMMARY: ", sep="")); print(length(which(duplicated(final_data$SYMBOL))))
		print("###############################################")
	
		#return
		return(final_data)
	}


	##apply the function
	#extract gene expression across studies and then save as a list preserving the labels
	list_final_datasets = llply(.data=bat_studies_considered, .fun=extract_expression_data)
		#llply is equivalent to lapply except that it will preserve labels and can display a progress bar.
	
	#see
	print("###############################################")
	print(paste("SEE DATASETS: ", paste(bat_studies_considered, collapse=" | "), sep="")); print(str(list_final_datasets))
	print("###############################################")

	
	##extract the names of all the genes for which we have expression data
	#extract the gene symbols of each data.frame in the list and save as a vector
	list_genes_expression_data = as.vector(unlist(sapply(X=list_final_datasets, "[", "SYMBOL")))
	
	#see
	print("###############################################")
	print(paste("CHECK NUMBER OF GENES ACROSS DATASETS: ", paste(bat_studies_considered, collapse=" | "), sep="")); print(length(list_genes_expression_data) == sum(sapply(X=list_final_datasets, nrow)))
	print("###############################################")

	#extract the unique cases
	list_genes_expression_data_unique = unique(list_genes_expression_data)
	print("###############################################")
	print(paste("UNIQUE GENES FOR WHICH WE HAVE EXPRESSION DATA: ", paste(bat_studies_considered, collapse=" | "), sep="")); print(length(list_genes_expression_data_unique))
	print("###############################################")


	
	#############################################################################
	####################### SELECT HIGHLY EXPRESSED GENES #######################
	#############################################################################
	
	#I think it is better to select those genes more expressed than compare expression between BAT and controls. In some studies the control is WAT of the same subject, but in others studies the controls are not differentiated cells... so I do not think it is a good idea to combine differential expression between studies. This is the cleanest way to combine multiple studies.
	
	#in case you need to do differential expression analyses, look:
		#ArrayExpress tutorial
			#https://www.bioconductor.org/packages/release/bioc/vignettes/ArrayExpress/inst/doc/ArrayExpress.pdf
		#very complete tutorial about array expression analyses
			#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#12_Linear_models
		#book chapter differential expression between pops
			#https://link.springer.com/content/pdf/10.1007%2F0-387-29362-0.pdf
	
	#it seems there are also option to merge array datasets
		#MergeMaid is an old package, but maybe there are novel options
	
	##create a function
	#for debugging: selected_dataset=list_final_datasets[[1]]; threshold_level=0.95
	high_expression = function(selected_dataset, threshold_level){
	
		#extract the expression value corresponding with the selected threshold
		selected_threshold = quantile(selected_dataset$average_gene_expression, probs=threshold_level)
	
		#select those rows of the selected dataset
		high_expression_subset = selected_dataset[which(selected_dataset$average_gene_expression > selected_threshold),]
	
		#check
		print("###############################################")
		print(paste("WE HAVE CORRECTLY SELECTED HIGHLY EXPRESSED GENES?", sep=""))
		print("###############################################")
		print(length(which(high_expression_subset$average_gene_expression<=selected_threshold)) == 0)
	
		#extract the gene symbols
		selected_genes = high_expression_subset$SYMBOL
	
		#return the results
		return(selected_genes)
	}
	
	
	##apply the function 
	#apply across the list of datasets and get a list as result
	genes_highly_expressed = llply(list_final_datasets, high_expression, threshold_level=0.95)
		#llply is equivalent to lapply except that it will preserve labels and can display a progress bar.
	
	#convert the list into a vector to have highly expressed genes across studies in the same vector
	genes_highly_expressed = as.vector(unlist(genes_highly_expressed))
	
	#select those unique cases
	genes_highly_expressed_unique = unique(genes_highly_expressed)
	
	#see
	print("###############################################")
	print(paste("NUMBER OF UNIQUE HIGHLY EXPRESSED GENES: ", paste(bat_studies_considered, collapse=" | "), sep="")); print(length(genes_highly_expressed_unique))
	print("###############################################")
	
	
	
	######################################################################
	####################### ANALYZE BAT CONNECTOME #######################
	######################################################################
	
	#load the information about BAT relationships
	bat_relationship = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/connectome_results/tables/appendix_S1_ordered.csv", sep=",", header=TRUE)
	#str(bat_relationship)
	#summary(bat_relationship)
	#head(bat_relationship) #we take the supplementary file
	
	#load the new file with the bat relationships obtained in 2020 from the original file of Jose, to check we are using the correct file
	bat_relationship_2 = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/connectome_results/tables/bat_relationship_check_2020.csv", sep=",", header=TRUE)
	#set as zero those cases with 2. I think these are cases not associated but that interact with other proteins that interact with UCP1
	bat_relationship_2[which(bat_relationship_2$BAT.relationship == 2),]$BAT.relationship <- 0
	#set the rows in alphabetic order according to gene name
	bat_relationship_2 = bat_relationship_2[order(bat_relationship_2$Genes),]
	#reset the row names
	row.names(bat_relationship_2) <- 1:nrow(bat_relationship_2)
	#add a new level to the Genes factors
	bat_relationship_2$Genes = factor(bat_relationship_2$Genes, c(unique(bat_relationship_2$Genes), "NRIP1"))
	#set "NRIP1 o RIP140" as NRIP1
	bat_relationship_2[which(bat_relationship_2$Genes == "NRIP1 o RIP140"),]$Genes <- "NRIP1"
	#remove the not used levels
	bat_relationship_2$Genes = droplevels(bat_relationship_2$Genes)
	
	#check
	print("###############################################")
	print(paste("WE HAVE CORRECT BAT CONNECTOME DATA ", paste(bat_studies_considered, collapse=" | "), sep="")); print(summary(as.vector(bat_relationship$Genes) == as.vector(bat_relationship_2$Genes))); print(summary(bat_relationship$BAT.relationship == bat_relationship_2$BAT.relationship))
	print("###############################################")
	
	
	#select genes based on the connectome status
	all_bat_genes = bat_relationship[which(bat_relationship$BAT.relationship %in% c(0,1)),]$Genes
	known_bat_genes = bat_relationship[which(bat_relationship$BAT.relationship %in% c(1)),]$Genes
	unknown_bat_genes = bat_relationship[which(bat_relationship$BAT.relationship %in% c(0)),]$Genes
	
	#see
	print("###############################################")
	print(paste("FOR HOW MANY BAT CONNECTOME GENES WE HAVE EXPRESSION DATA: ", paste(bat_studies_considered, collapse=" | "), sep="")); print(paste(length(which(all_bat_genes %in% list_genes_expression_data_unique)), " out of ", length(all_bat_genes), sep=""))
	print("###############################################")

	
	## select all genes outside the connectome (both known and unknown to be associated with BAT)
	#load the connectome with UCP1 as core gene
	ucp1_conn = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/data/human_connectome/UCP1.txt", sep="\t", header=T)
	#str(ucp1_conn)
	#head(ucp1_conn)
	#summary(ucp1_conn)
	
	#create a vector with all coding human genes ordered alphabetically from the UCP1 connectome
	all_genes = sort(ucp1_conn$Target)
	
	#select all genes except those included in the connectome
	random_genes = all_genes[which(!all_genes %in% all_bat_genes)]
		#it is not specially relevant if we remove all the BAT connectome genes or just the known BAT genes for the test of known, and just unknown BAT genes for the test of unknown. We are talking of 168 genes among 16,000. We can just remove all the BAT connectome genes (as we did for the previous BAT-known analyses) without impacting significantly the pool of control genes.
		#in the same vein, yes, the same gene could be considered twice as control, but the probability is very low given that we are selecting around 100 genes each time within a pool of 16000 genes.
	
	#see
	print("###############################################")
	print(paste("FOR HOW MANY RANDOM GENES WE HAVE EXPRESSION DATA: ", paste(bat_studies_considered, collapse=" | "), sep="")); print(paste(length(which(random_genes %in% list_genes_expression_data_unique)), " out of ", length(random_genes), sep=""))
	print("###############################################")


	##calculate the number of highly expressed BAT genes in each set of the BAT connectome
	number_all_bat_genes = length(which(all_bat_genes %in% genes_highly_expressed_unique))
	number_known_bat_genes = length(which(known_bat_genes %in% genes_highly_expressed_unique))
	number_unknown_bat_genes = length(which(unknown_bat_genes %in% genes_highly_expressed_unique))
	
	
	##calculate the number of highly expressed BAT genes in random sets
	#set the seed
	set.seed(98743)
	
	#open empty vectors
	number_random_all_bat_genes = NULL
	number_random_known_bat_genes = NULL
	number_random_unknown_bat_genes= NULL
	
	#calculate 1 million random sets of genes
	for(i in 1:100000){
		
		#select a number of random genes matching each BAT connectome set
		random_all_bat_genes = sample(1:length(random_genes), length(all_bat_genes), replace=FALSE)
		random_known_bat_genes = sample(1:length(random_genes), length(known_bat_genes), replace=FALSE)
		random_unknown_bat_genes = sample(1:length(random_genes), length(unknown_bat_genes), replace=FALSE)
			#We are using no replacement, because we want every time a gene is included in the random set, then that gene cannot be selected again. In other words, we do not want the same gene two times in the same control set. This is unlikely given the great pool of controls, but we use this just in case.
	
		#extract the number of random genes in the highly expressed set
		number_random_all_bat_genes = append(number_random_all_bat_genes, length(which(random_genes[random_all_bat_genes] %in% genes_highly_expressed_unique))) 
		number_random_known_bat_genes = append(number_random_known_bat_genes, length(which(random_genes[random_known_bat_genes] %in% genes_highly_expressed_unique))) 
		number_random_unknown_bat_genes = append(number_random_unknown_bat_genes, length(which(random_genes[random_unknown_bat_genes] %in% genes_highly_expressed_unique))) 
	}
	
	
	## calculate the pvalues
	
	#Calculate p.vale as the proability of a random set of genes will have more BAT highly expressed genes than the number found in the BAT connectome
	pval_all_bat_genes = prop.test(x=length(which(number_random_all_bat_genes >= number_all_bat_genes)), n=length(number_random_all_bat_genes)) 
		#https://stats.stackexchange.com/questions/167164/test-of-equal-proportions-with-zero-successes
	
	#Calculate p.vale as the proability of a random set of genes will have more BAT highly expressed genes than the number found in the known bat genes
	pval_known_bat_genes = prop.test(x=length(which(number_random_known_bat_genes >= number_known_bat_genes)), n=length(number_random_known_bat_genes)) 
		#https://stats.stackexchange.com/questions/167164/test-of-equal-proportions-with-zero-successes
	
	#Calculate p.vale as the proability of a random set of genes will have more BAT highly expressed genes than the number found in the unknown bat genes
	pval_unknown_bat_genes = prop.test(x=length(which(number_random_unknown_bat_genes >= number_unknown_bat_genes)), n=length(number_random_unknown_bat_genes)) 
		#https://stats.stackexchange.com/questions/167164/test-of-equal-proportions-with-zero-successes
	
	#print
	print("###############################################")
	print(paste("P-VALUE BAT CONNECTOME: ", paste(bat_studies_considered, collapse=" | "), sep=""))
	print("###############################################")
	print(pval_all_bat_genes)
	print("###############################################")
	print(paste("P-VALUE KNOWN BAT GENES: ", paste(bat_studies_considered, collapse=" | "), sep=""))
	print("###############################################")
	print(pval_known_bat_genes)
	print("###############################################")
	print(paste("P-VALUE UNKNOWN BAT GENES: ", paste(bat_studies_considered, collapse=" | "), sep=""))
	print("###############################################")
	print(pval_unknown_bat_genes)

	#return the names of the datasets processed
	return(bat_studies_considered)
}



#####################################
###### PARALLELIZE THE PROCESS ######
#####################################

#BAT studies
sets_of_datasets = list(ids_bat_studies[which(ids_bat_studies %in% c("E-GEOD-27657", "E-GEOD-54280"))], ids_bat_studies)

#set up cluster
clust <- makeCluster(length(sets_of_datasets), outfile="") #outfile let you to see the output in the terminal "https://blog.revolutionanalytics.com/2015/02/monitoring-progress-of-a-foreach-parallel-job.html"
registerDoParallel(clust)

#run the function
sets_of_datasets_processed = foreach(i=sets_of_datasets, .packages=c("ArrayExpress", "arrayQualityMetrics", "dplyr", "plyr", "hugene10sttranscriptcluster.db", "hgu133plus2.db")) %dopar% {
    full_analysis(bat_studies_considered=i)
}

#stop the cluster 
stopCluster(clust)

#check
print("###############################################")
print(paste("WE PROCESSED ALL SELECTED DATASETS?", sep="")); print(identical(sets_of_datasets, sets_of_datasets_processed))
print("###############################################")
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



########################################################################
####################### CALCULATE GEN EXPRESSION #######################
########################################################################


##write a function to do that
#selected_ids_bat_studies=ids_bat_studies[1] #for debugging
extract_expression_data = function(selected_ids_bat_studies){

	#make a directory for the selected dataset
	system(paste("mkdir -p array_express_gene_expression/", selected_ids_bat_studies, sep=""))
		#use mkdir -p because if you don't and run again, you will get an error that it exists. "p" flags is for "no error if existing, make parent directories as needed"

	#see one specific dataset
	print("###############################################")
	print(paste("STARTING", selected_ids_bat_studies, ": ", sep=""))
	print("###############################################")


	#load the AE set filtered
	AEsetnorm_filter = readRDS(paste("array_express_raw/", selected_ids_bat_studies, "/", selected_ids_bat_studies, "_expression_filter_1.Rds", sep=""))


	##prepare gene expression data
	#get a data.frame with rows representing genes and columns representing arrays. To know which columns go with that tissue or experimental treatment, we can rely on the phenoData information inherited from AEset_human.
	expression_matrix = data.frame(exprs(AEsetnorm_filter))
		#access the expression and error measurements of assay data stored in an object derived from the ‘eSet-class’.
		#head(expression_matrix)

	#we have a value of expression per sample, so we have to find a way to summarize and get a value across all samples per gene. We have filtered the sample to be considered in a previous script, so we can just use the median to summarize the gene expression across all the remaining samples. These should be only BAT-related samples of humans thanks to the filtering. In the same experiments we can have BAT sample of males and females of different ages, or brown adipocytes generated in different ways. But thanks to the filtering, we know that the remaining samples belong to experimental groups that show BAT-like features.
	#We will use median like in the rest of the paper in order to summarize.
	
	#calculate the median expression across all samples per prob
	median_expression = as.data.frame(apply(X=expression_matrix, MARGIN=1, FUN=median))
	#set the column name as average expression
	colnames(median_expression) = "average_gene_expression"
	str(median_expression)
	#check
	print("###############################################")
	print(paste("MERGIN", selected_ids_bat_studies, " OK?: ", sep="")); print(identical(row.names(median_expression), row.names(expression_matrix)))
	print("###############################################")

	#merge the result with the original data.frame using the row names 
	expression_matrix = merge(expression_matrix, median_expression, by="row.names")

	#select only the row names and the average gene expression
	expression_matrix_average = expression_matrix[,c("Row.names", "average_gene_expression")]
	#str(expression_matrix_average)


	##gene annotation
	
	#We used the function select from AnnotationDbi to query the gene symbols and associated short descriptions for the transcript clusters. For each cluster, we added the gene symbol (SYMBOL) and a short description of the gene the cluster represents (GENENAME).
		#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#11_Annotation_of_the_transcript_clusters	

	#select the AnnotationDb object
	#we have to convert the probs ID to the gene symbols used by Yuval. In all cases, we select the annotation package that it is not probset. I understand that we do not need probset summarization, which leads to exons, but gene summarization
		#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#74_Old_and_new_%E2%80%9Cprobesets%E2%80%9D_of_Affymetrix_microarrays
	if(AEsetnorm_filter@annotation == "pd.hugene.1.0.st.v1"){

		#we select the hugene 10 annotation package
		require(hugene10sttranscriptcluster.db)
		selected_annotation = hugene10sttranscriptcluster.db
			#http://bioconductor.org/packages/release/BiocViews.html#___AnnotationData
			#http://bioconductor.org/packages/release/data/annotation/manuals/hugene10sttranscriptcluster.db/man/hugene10sttranscriptcluster.db.pdf
	}
	if(AEsetnorm_filter@annotation == "pd.hg.u133.plus.2"){
		
		#we select the hg 133 annotation package
		require(hgu133plus2.db)
		selected_annotation = hgu133plus2.db
			#http://bioconductor.org/packages/release/BiocViews.html#___AnnotationData
			#http://bioconductor.org/packages/release/data/annotation/manuals/hgu133plus2.db/man/hgu133plus2.db.pdf
	}

	#using the selected annotation package, extract the gene symbols for the probs we have (rows in the expression data)
	gene_symbols_probs = mapIds(x=selected_annotation, keys=expression_matrix_average$Row.names, column=c("SYMBOL"), keytype="PROBEID", multiVals="first")
		#multivals:
			#first: This value means that when there are multiple matches only the 1st thing that comes back will be returned. This is the default behavior. 


		#https://support.bioconductor.org/p/69378/#69379
		#https://support.bioconductor.org/p/70769/


	#anno_palmieri <- AnnotationDbi::select(hugene10sttranscriptcluster.db, keys = (featureNames(AEsetnorm_filter)), columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")
	

	#anno_palmieri <- subset(anno_palmieri, !is.na(SYMBOL))


	#Traditionally, Affymetrix arrays (the so-called 3’ IVT arrays) were probeset based: a certain fixed group of probes were part of a probeset which represented a certain gene or transcript (note however, that a gene can be represented by multiple probesets). The more recent “Gene” and “Exon” Affymetrix arrays are exon based and hence there are two levels of summarization to get to the gene level. The “probeset” summarization leads to the exon level. The gene / transcript level is given by “transcript clusters”. Hence, the appropriate annotation package for our chip type is called hugene10sttranscriptcluster.db.
	#On the left side, we see plenty of probes for each Exon / probeset (i.e. each colour): therefore, a summarization on the probeset / exon level makes sense. In the gene type array, however, only a small proportion of the original probes per probeset is included. Thus, a summarization on the probeset / exon level is not recommended for “Gene” arrays but nonetheless possible by using the hugene10stprobeset.db annotation package. Note that furthermore, there are also no longer designated match/mismatch probes present on “Gene” and “Exon” type chips. The mismatch probe was initially intended as base-level for background correction, but hasn’t prevailed due to more elaborate background correction techniques that do not require a mismatch probe.
		#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#74_Old_and_new_%E2%80%9Cprobesets%E2%80%9D_of_Affymetrix_microarrays

	#HAY QUE REVISAR QUE EL MISMO TRANSCRITO NO ESTÁ ASOCIADO A DIFERENTES GENES, ESOS CASOS HAY QUE QUITARLOS
		#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#111_Removing_multiple_mappings








gene_symbols_probs = gene_symbols_probs[which(!is.na(gene_symbols_probs))] #CHECK HOW MANY ARE LOST
gene_symbols_probs = as.data.frame(gene_symbols_probs)
gene_symbols_probs$PROBEID = row.names(gene_symbols_probs)
str(gene_symbols_probs)


require(dplyr)
anno_grouped <- group_by(gene_symbols_probs, PROBEID)
anno_summarized <- dplyr::summarize(anno_grouped, no_of_matches = n_distinct(gene_symbols_probs))
head(anno_summarized)
summary(anno_summarized)


merged_data = merge(expression_matrix_average, gene_symbols_probs, by="row.names") #CHECK THIS
str(merged_data)
head(merged_data)
summary(merged_data)


merged_data[which(duplicated(merged_data$gene_symbol)),]


merged_data_aggregated = aggregate(average_gene_expression ~ gene_symbols_probs, merged_data, median)
	#some recommend to summarize probs by gene using mean or median. We have used median to summarize across the manuscript
		#https://www.biostars.org/p/271379/#271588
		#https://www.biostars.org/p/336130/

	#OJO
		#If they are different isoforms, doing an average of two might not be appropriate cause one might no be expressed and that brings gene expression down..
		#Yes, of course, in which case you will have to write some code to check for these situations in which the expression is so low for one probe such that it is negligible.
			#https://www.biostars.org/p/271379/#271588

		#otros dicen the usar el valor maximo
			#https://support.bioconductor.org/p/70133/


merged_data_aggregated[which(duplicated(merged_data_aggregated$gene_symbol)),]


merged_data_aggregated_ordered = merged_data_aggregated[order(merged_data_aggregated$average_gene_expression, decreasing=TRUE),]

merged_data_aggregated_high_expression = merged_data_aggregated[which(merged_data_aggregated$average_gene_expression > quantile(merged_data_aggregated$average_gene_expression, probs=0.75)),]

#I think it is better to select those genes more expressed than compare with controls. In some studies the control is WAT of the same subject, but in others studies the control are not differentiated cells... so I do not think it is a good idea to combine differential expression between studies. This is the cleanest way to combine multiple studies.




	#MIRA
	#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#12_Linear_models


	#book chapter differential expression between pops
		#https://link.springer.com/content/pdf/10.1007%2F0-387-29362-0.pdf


#ESTAS LINEAS SON DEL ARRAYEXPRESS TUTORIAL
#Now that we have ensured that the data are well processed, we can search for differentially expressed genes using the package limma. To understand the details of each steps, please see the limma user guide.
library("limma")
facs = pData(AEsetnorm_filter)[,column_filter_names]
f = factor(facs)
design = model.matrix(~0+f)
colnames(design) = levels(f)
fit = lmFit(AEsetnorm_filter, design)
cont.matrix = makeContrasts(BATvsWAT = BAT-WAT, levels=design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

#Here we end up with a list of genes that are differentially expressed between BAT and WAR, but one can perform other comparisons with the same data.
res = topTable(fit2, coef = "BATvsWAT", adjust = "BH")

#This could now be followed by an integrative analysis of the data, a complex and open-ended task for which essential tools are provided in the Bioconductor project: the quality of the datasets could be assessed with the help of the arrayQualityMetrics package (Kauffmann et al., 2009), they could be normalized and analysed for differential expression of genes and gene sets (Hahne et al., 2008), and the combination of different datasets is facilitated, for example, by the MergeMaid package (Cope et al., 2004).






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


#load the connectome with UCP1 as core gene
ucp1_conn = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/data/human_connectome/UCP1.txt", sep="\t", header=T)
str(ucp1_conn)
head(ucp1_conn)
summary(ucp1_conn)


## select all genes outside the connectome (both known and unknown to be associated with BAT)
#create a vector with all codifican human genes ordered alphabetically from the UCP1 connectome
all_genes = sort(ucp1_conn$Target)
#select all genes except those included in the connectome
random_genes = all_genes[-which(all_genes %in% bat_relationship[which(bat_relationship$BAT.relationship %in% c(0,1)),]$Genes) ]


#threshold_top = 5000 #we cannot compare the median expression between sets, because we are going to have mutliple studies and I think it is not ok to combine arrays of different studies. Remember, that in each study, we compare the intensity signal of each array in order to detect outliers, e.g., arrays with a different level of background noise, which could bias results. we should do that with the arrays of all studies in order to combine them.


#VAMOS A MIRAR DIRECTMANETE CUANTOS DE LOS CANDIDATOS SE EXPRESAN EN EL BAT, HAY QUE BUSCAR UN CUTOFF Y ELIMINAR LOW EXPRESSED GENES, RPIMER CUARTIL DE MAS EXPRESADOS?. El revisor 2 no dice nada de expression diferencial, solo que analizemos transcriptoma data y validemos los candidatos como BAT markers. Necesitamos mostrar que hay candidatos altamente expresados en el BAT.
#perdemos poder estadistico, pero ganamos mas finura en la parte biologica al mirar especificamente WAT vs BAT. El conectoma está vlaidadio estadisticamente, vmaos a hora a añadir info biologica.
#si sale p-value, bien, si no, decimos que tenemos unos cuantos con eivdencia de upregulation y que aun asi los que no tenemos eviedncia pueden serlo, mira thyroid hormones. simplememte refinamos aun mas la lista.


#number_genes_candidate = length(which(unknown_bat_genes %in% merged_data_aggregated_ordered[1:threshold_top, "gene_symbols_probs"]))
number_genes_candidate = length(which(unknown_bat_genes %in% merged_data_aggregated_high_expression$gene_symbols_probs))


#NO VAMOS A HACER TEST, PORQUE VA CAMBIANDO SEGUN LOS ESTUDIOS QUE INCLUIMOS, HAY QUE TENER EN CUENTA QUE NO TENEMOS UNA CURACION DE TODOS LOS ESTUDIOS ANALIZANDO BAT... ESTE TEJIDO NO ESTA EN GTEEX... ASI QUE MEJOR DECIR QUE TENEMOS XX CANDIDATOS ENTRE EL PRIMER CUARTIL  DE TODOS LOS ESTUDIOS CONSIDERADOS Y PUNTO... ASI SI NOS PIDEN MAS ESTUDIOS, ESO NO VA A CAMBIAR...

number_genes_random = NULL
for(i in 1:100000){
	random_expression_rows = sample(1:length(random_genes), length(unknown_bat_genes), replace=FALSE)
	#number_genes_random = append(number_genes_random, length(which(random_genes[random_expression_rows] %in% merged_data_aggregated_ordered[1:threshold_top, "gene_symbols_probs"]))) #CHECK THE NUMBER OF RANDOM GENES WITH EXPRESSION DATA
	number_genes_random = append(number_genes_random, length(which(random_genes[random_expression_rows] %in% merged_data_aggregated_high_expression$gene_symbols_probs))) #CHECK THE NUMBER OF RANDOM GENES WITH EXPRESSION DATA
}


#Calculate p.vale as the proability of a random pair of genes will have lower biological distance than the median of biological distance between all UCP1 connectome genes.
pval = prop.test(x=length(which(number_genes_random >= number_genes_candidate)), n=length(number_genes_random)) #https://stats.stackexchange.com/questions/167164/test-of-equal-proportions-with-zero-successes
pval









#THE FINAL OUTPUT OF THIS SCRIPT IS
	#a tab delimited file for each expression dataset, having two columns, the gene symbol and the average gene expression across subjects, treatments and probs. BUT ONLY FOR BAT NOT WAT.
	#the rows will be in decreasing order following the expression levels


#IN THE NEXT SCRIPT
#load a list with all gene names included in the human genome connectome and those that are BAT candidates genes.
#create a function that take the top as an argument, for example top 100 of genes more expressed
	#create a function that
		#read only the 100 first rows of each dataset, which are ordered in decreasing order based on gene expression
		#vind all rows indicating the name of the dataset
	#in the final file, select the unique gene symbols. These are the genes present in the top 100 of any of the studies
	#create a loop with 10,001 iterations
		#in the first iteration select the BAT candidates and calculate how many of them are included in the list of top genes previously created. Save the number.
		#for the rest iterations, create a random set of gene names with the same size than BAT candidates we have. In each case, calculate how many are in the top across studies, and then save the number.
	#calculate a p-value comparing the number of BAT candidates in the expression top with that of random genes.
#this can done for 2000 to 10 tops and calculate and enrichment curve

}
#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



#####################################################################################
########################## OBTAINING BAT GENE EXPRESSSION ###########################
#####################################################################################

#Script for calculating cleaned dataset with gene expression values in BAT

#Part of the code of this script comes from a tutorial of ArrayExpress
	#http://www.bioconductor.org/packages/release/bioc/vignettes/ArrayExpress/inst/doc/ArrayExpress.pdf



#################################################################
####################### REQUIRED PACKAGES #######################
#################################################################

require(ArrayExpress) #for loading data from ArrayExpress
require(arrayQualityMetrics) #for normalizing expression data in each dataset



#######################################
#### CONNECTING WITH ARRAY EXPRESS ####
#######################################

##connect with ArrayExpress
#ArrayExpress is a public repository for transcriptomics and related data. Among other data includes micrarrays, which is a substrate to which are attached millions of single-stranded (c)DNA complementary to the (c)DNA you wish to detect
	#https://online.stat.psu.edu/stat555/node/28/
	#https://www.ebi.ac.uk/training/online/courses/array-express-discover-functional-genomics-data-quickly-and-easily/#vf-tabs__section--contents

#set working directory to results, because ArrayExpress save stuff during the analyses automatically
setwd("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/results_2022/array_express_raw")

#make a query for BAT
sets=queryAE(keywords="brown+adipose+tissue", species="homo+sapiens")
	#keywords: the keyword(s) of interest. To use several words, they must be separated by a "+" as shown in the examples.
	#species: the specie(s) of interest.

#take a look
str(sets)
	#we get the same, 15 experiments, just like we search in the web, the same IDs
	#https://www.ebi.ac.uk/arrayexpress/search.html?query=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22


##select sets with raw data available
#I prefer to select only those with raw data so I can process the data in the same way for all datasets.
sets_raw=sets[which(sets$Raw=="yes"),]

#check
print("#########################################")
print("CHECK THAT ALL DATASETS HAVE DATA FOR HUMANS AND INCLUDE RAW DATA")
print(length(which(sets_raw$Raw != "yes"))==0 & length(which(!grepl("Homo sapiens", sets_raw$Species)))==0)
print("#########################################")


##select only those experiments measuring gene expression at BAT
#sets with BAT data
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8564/?keywords=brown+adipose+tissue&organism=Human&exptype%5B%5D=&exptype%5B%5D=&array=
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-4031/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-56635/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2602/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-57896/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-49795/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-27657/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-19643/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=

#brown preadipocytes
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-68544/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=

#it seems but not sure
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-54280/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-25/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=

#strange case report of BAT in visceral fat
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-49795/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=

#set with expression done on WAT before and after cold exposure
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-67297/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=

#create a vector with the IDs
ids_bat_studies = c("E-MTAB-25", "E-GEOD-56635")
	#E-MTAB-25 fail, Error in `.rowNamesDF<-`(x, value = value) :  duplicate 'row.names' are not allowed
		#maybe: https://www.biostars.org/p/62988/


#######################################
#### OBTAIN DATA FOR EACH DATASET #####
#######################################

#extend the time window during which a function can download, for example, download.file, which is used by ArrayExpress function
getOption('timeout')
options(timeout=1000)
getOption('timeout')
	#https://stackoverflow.com/questions/35282928/how-do-i-set-a-timeout-for-utilsdownload-file-in-r


##write a function to do that
#selected_ids_bat_studies=ids_bat_studies[1] #for debugging
extract_expression_data = function(selected_ids_bat_studies){

	#see one specific dataset
	print("###############################################")
	print(paste("INFO", selected_ids_bat_studies, ": ", sep=""))
	print(sets_raw[which(sets_raw$ID == selected_ids_bat_studies), 1:7])
	print("###############################################")

	#Use ArrayExpress for loading the selected dataset
	selected_aeset=ArrayExpress(selected_ids_bat_studies)
		#If the identifier refers to an Affymetrix experiment, the output is an AffyBatch, if it refers to a one-colour experiment using a platform other than Affymetrix, the output is an ExpressionSet. The ArrayExpress function extracts feature intensity summaries from columns of the raw data files based on the common conventions for the data file sources. If the data source is not recognized, or the file does not have the expected column names, the user is asked to explicitly provide the name of the column(s) to extract, for instance, ‘Cy3 Median’. In some cases, there is a mismatch between the sample or feature annotations and the intensity data files; in such cases, a warning is emitted, the phenoData and/or featureData components are left empty and an incomplete (but syntactically valid) object is returned. 

	#if the dataset is a list because we have data from different species
	if(is.list(selected_aeset) & length(selected_aeset)>1){

		#create a function to extract the annotation of each dataset
		get_annotation <- function(my_list){return(my_list@annotation)}
			#https://stackoverflow.com/questions/24611365/subsetting-slot-from-a-list-of-data-frames

		#get the annotation
		annotations_list = lapply(selected_aeset, get_annotation)

		#extract the name of the dataset that has human genes as annotation
		name_human_dataset = names(annotations_list[which(grepl("hugene", annotations_list, fixed=TRUE))])

		#select the corresponding dataset
		AEset_human = selected_aeset[name_human_dataset]

		#only one element should be selected so we select it and remove the list envelop
		AEset_human = AEset_human[[1]]
	} else { #if not
		
		#do nothing
		AEset_human = selected_aeset
	}
	
	#make some checks to avoid datasets without human data
	if(!grepl("hugene", AEset_human@annotation, fixed=TRUE)){
		print("###############################################")
		stop(paste(selected_ids_bat_studies, " HAS NO HUMAN GENES", sep=""))
		print("###############################################")
	}

	#make some checks to avoid datasets that are not affymetrix experiments
	if(AEset_human@manufacturer != "Affymetrix"){
		print("###############################################")
		stop(paste(selected_ids_bat_studies, " IS NOT AN AFFYMETRIX EXPERIMENT", sep=""))
		print("###############################################")
	}

	#make some checks to avoid datasets that are not gene or expression sets
	if(class(AEset_human) != "GeneFeatureSet" & class(AEset_human) != "ExpressionFeatureSet"){
		print("###############################################")
		stop(paste(selected_ids_bat_studies, " IS NOT CLASS GENE OR EXPRESSION FEATURE SET", sep=""))
		print("###############################################")
	}

	#According the ArrayExpress vignette, we can use RMA normalisation to process the raw data, please read the rma function help.
	AEsetnorm = oligo::rma(AEset_human, background=TRUE, normalize=TRUE, target="core")
		#the tutorial says affy:rma, but for me it does not work. I have to use oligo::rma.
			#https://support.bioconductor.org/p/71097/
		#oligo::rma
			#Robust Multichip Average preprocessing methodology. This strategy allows background subtraction, quantile normalization and summarization (via median-polish). In other words, during the RMA normalisation, a Tukey's 'median polish' is applied using information from all probesets.
			#object: Exon/HTA/Expression/Gene/SnpCnv-FeatureSet object.
			#background: Logical - perform RMA background correction?
			#normalize: Logical - perform quantile normalization?
			#target: Level of summarization (e.g., at the gene level... only for Exon/Gene arrays)
				#core would summarize to the gene level. However, if your array is an 'Exon' array, this may still only summarise to exon-level, requiring you to perform a further [manual] summarisation to gene-level.
					#https://www.biostars.org/p/271379/#271588


#see the data 
head(pData(AEsetnorm))
	#this is the same found in detailed sample information
		#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-67297/samples/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=

#To check the normalisation efficiency, we can run a quality assessment. For details on the arguments used, please read the arrayQualityMetrics vignette.
#get the columns with the experimental factors (e.g., age, sex, treatment, cell type..)
fac = grep("Sample_source|sex", colnames(pData(AEsetnorm)), value=T)
	#we HAVE FACTOR.VALUE CELL TYPE OR ORGANISM
qanorm = arrayQualityMetrics(AEsetnorm, outdir="QAnorm", intgroup=fac, force=TRUE)

#plot
#pdf("/home/dftortosa/Desktop/eso.pdf")
#qanorm
#dev.off()

#Now that we have ensured that the data are well processed, we can search for differentially expressed genes using the package limma. To understand the details of each steps, please see the limma user guide.
library("limma")
facs = pData(AEsetnorm)[,fac]
facs[which(grepl("neck", facs[,1])), 1]="N"
facs[facs[,1]!="N", 1]="S"
facs[facs[,2]=="male", 2]="male"
facs[facs[,2]=="female", 2]="female"
facs = paste(facs[,2], facs[,1], sep=".")
f = factor(facs)
design = model.matrix(~0+f)
colnames(design) = levels(f)
fit = lmFit(AEsetnorm, design)
cont.matrix = makeContrasts(male.NvsS = male.N-male.S, female.NvsS = female.N-female.S, Diff=(female.N-female.S)-(male.N-male.S), levels=design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

#Here we end up with a list of genes that are differentially expressed between neck and subcutaneous adipose tissue in males and females, respectively, but one can perform other comparisons with the same data. For example, differences across age
res_male = topTable(fit2, coef = "male.NvsS", adjust = "BH")
res_female = topTable(fit2, coef = "female.NvsS", adjust = "BH")

#This could now be followed by an integrative analysis of the data, a complex and open-ended task for which essential tools are provided in the Bioconductor project: the quality of the datasets could be assessed with the help of the arrayQualityMetrics package (Kauffmann et al., 2009), they could be normalized and analysed for differential expression of genes and gene sets (Hahne et al., 2008), and the combination of different datasets is facilitated, for example, by the MergeMaid package (Cope et al., 2004).


expression_matrix = data.frame(exprs(AEsetnorm))
head(expression_matrix)

#we have a value of expression per subject, so we have to find a way to summarize and get a value across all subject per gene. We would use median like in the rest of the paper
	#it makes sense to average expression across sex? but what about case/controls? or cold vs non cold?

expression_matrix$average_gene = data.frame(apply(X=expression_matrix, MARGIN=1, FUN=median))
str(expression_matrix)

expression_matrix_average = expression_matrix[,"average_gene"]
colnames(expression_matrix_average) = "average_gene_expression"

#we have to convert the probs ID of affy_hugene_1_0_st_v1 to the gene symbols used by Yuval
	#https://www.biostars.org/p/69597/
	#https://support.bioconductor.org/p/42839/

require(hugene10sttranscriptcluster.db) #I guess hu gene 10 is hugene_1_0... I ahve seen the same for hugene20 as hugene_2_0
tail(keys(hugene10sttranscriptcluster.db, keytype="SYMBOL"))
k <- keys(hugene10sttranscriptcluster.db, keytype="PROBEID")
select(hugene10sttranscriptcluster.db, keys=k, columns=c("SYMBOL","GENENAME"), keytype="PROBEID")

gene_symbols_probs = mapIds(hugene10sttranscriptcluster.db, keys=k, column=c("SYMBOL"), keytype="PROBEID", multiVals="first")
	#mapIds deals with duplicates...
		#https://support.bioconductor.org/p/69378/#69379

		#https://support.bioconductor.org/p/70769/

gene_symbols_probs = gene_symbols_probs[which(!is.na(gene_symbols_probs))]
gene_symbols_probs = as.data.frame(gene_symbols_probs)
str(gene_symbols_probs)

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


merged_data_aggregated[order(merged_data_aggregated$average_gene_expression, decreasing=TRUE),]


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


merged_data_aggregated[which(merged_data_aggregated$gene_symbol %in% unknown_bat_genes),]








require(biomaRt)
# replace the affyID with gene symbol
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="http://uswest.ensembl.org/", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

hgnc <- getBM(attributes = c("affy_hugene_1_0_st_v1", "hgnc_symbol","ensembl_gene_id","entrezgene","chromosome_name","start_position","end_position","band"), filters = "affy_hugene_1_0_st_v1", values=tab$ID, mart = mart)
	#pd.hugene.1.0.st.v1

# Now match the array data probesets with the genes data frame
m <- match(as.numeric(tab$ID), hgnc$affy_hugene_1_0_st_v1)
# And append e.g. the HGNC symbol to the array data frame
tab$hgnc <- hgnc[m, "hgnc_symbol"]
	#https://www.biostars.org/p/69597/


#there is an annotation tool for ugene_1_0_st_v1
	#http://www.bioconductor.org/packages/release/data/annotation/html/pd.hugene.1.0.st.v1.html

data(pmSequence) #make a try


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
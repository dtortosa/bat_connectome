#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



#####################################################################################
########################## OBTAINING BAT GENE EXPRESSSION ###########################
#####################################################################################

#Script for calculating cleaned dataset with gene expression values in BAT

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



#######################################
#### CONNECTING WITH ARRAY EXPRESS ####
#######################################

##connect with ArrayExpress
#ArrayExpress is a public repository for transcriptomics and related data. Among other data includes micrarrays, which is a substrate to which are attached millions of single-stranded (c)DNA complementary to the (c)DNA you wish to detect
	#https://online.stat.psu.edu/stat555/node/28/
	#https://www.ebi.ac.uk/training/online/courses/array-express-discover-functional-genomics-data-quickly-and-easily/#vf-tabs__section--contents

#set working directory to results, because ArrayExpress save stuff during the analyses automatically
wd="/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/results_2022"
setwd(wd)

#create a directory to save the results and move to it
system(paste("mkdir -p array_express_raw", sep=""))
	#use mkdir -p because if you don't and run again, you will get an error that it exists. "p" flags is for "no error if existing, make parent directories as needed"
setwd("array_express_raw")

#make a query for BAT
sets=queryAE(keywords="brown+adipose+tissue", species="homo+sapiens")
	#keywords: the keyword(s) of interest. To use several words, they must be separated by a "+" as shown in the examples.
	#species: the specie(s) of interest.

#take a look
str(sets)
	#we get the same, 15 experiments, just like we search in the web, the same IDs
	#https://www.ebi.ac.uk/arrayexpress/search.html?query=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22


##select sets with raw data available
#I prefer to select only those with raw data so I can process the data in the same way for all datasets. For example, if a dataset has been processed using RMA, expression has been log transformed, but not if the MAS5 processing was applied. Therefore gene expression of each study is not comparable. In our case, we are going to calculate tops in each dataset, but I think it is better to maintain the same treatment just in case...
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
#selected_ids_bat_studies=ids_bat_studies[2] #for debugging
extract_expression_data = function(selected_ids_bat_studies){

	#make a directory for the selected dataset
	system(paste("mkdir -p ", selected_ids_bat_studies, sep=""))
		#use mkdir -p because if you don't and run again, you will get an error that it exists. "p" flags is for "no error if existing, make parent directories as needed"

	#see one specific dataset
	print("###############################################")
	print(paste("INFO", selected_ids_bat_studies, ": ", sep=""))
	print(sets_raw[which(sets_raw$ID == selected_ids_bat_studies), 1:7])
	print("###############################################")

	#Use ArrayExpress for loading the selected dataset
	selected_aeset=ArrayExpress(selected_ids_bat_studies, path=selected_ids_bat_studies, save=FALSE)
		#If the identifier refers to an Affymetrix experiment, the output is an AffyBatch, if it refers to a one-colour experiment using a platform other than Affymetrix, the output is an ExpressionSet. The ArrayExpress function extracts feature intensity summaries from columns of the raw data files based on the common conventions for the data file sources. If the data source is not recognized, or the file does not have the expected column names, the user is asked to explicitly provide the name of the column(s) to extract, for instance, ‘Cy3 Median’. In some cases, there is a mismatch between the sample or feature annotations and the intensity data files; in such cases, a warning is emitted, the phenoData and/or featureData components are left empty and an incomplete (but syntactically valid) object is returned.
		#arguments
			#accession: an ArrayExpress experiment identifier.
			#path: the name of the directory in which the files downloaded on the ArrayExpress repository will be extracted. The default is the current directory.
			#save: if TRUE, the files downloaded from the database will not be deleted from path after executing the function.
				#we select FALSE because we are just going to save a rds file with expression data after preprocessing
			#dataCols: by default, for the raw data, the columns are automatically selected according to the scanner type. If the scanner is unknown or if the user wants to use different columns than the default, the argument 'dataCols' can be set. For two colour arrays it must be a list with the fields 'R', 'G', 'Rb' and 'Gb' giving the column names to be used for red and green foreground and background. For one colour arrays, it must be a character string with the column name to be used. These column names must correspond to existing column names of the expression files.
				#I guess this is related with the columns where the colour of the array is presented. Remember that the result of array experiments is a color in different probs
					#https://online.stat.psu.edu/stat555/node/28/

	#load the data
	#I am saving the data, but right now I am not sure how to load the data to my R session.

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
	
	#make some checks to avoid datasets without HuGenes
	if(!grepl("hugene", AEset_human@annotation, fixed=TRUE)){
		print("###############################################")
		stop(paste(selected_ids_bat_studies, " HAS NO HUMAN GENES", sep=""))
		print("###############################################")
	} #According to some, if you have HuGene arrays, you should use RMA (or GCRMA) normalization, since you don't have mismatch probes (although you have other types of control probes).
		#https://www.biostars.org/p/211389/

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


	##RLE
	#Before doing the RMA processing, you can calculate for each array how many transcripts (probes) are deviated from the median across arrays (i.e., median expression of a probe across all arrays). You can plot it as a boxplot with boxes having a larger extension therefore indicate an unusually high deviation from the median in a lot of transcripts, suggesting that these arrays are different from most of the others in some way.
	#You can keep these samples in mind for heatmap cluster analysis later on in the workflow. Arrays that are confirmed to be outliers by heatmap analysis could be removed for subsequent analysis.
	#We are not going to do this, because we are already detecting arrays that are outliers with different plots and statistical tests using arrayQualityMetrics, see below.
		#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#8_Relative_Log_Expression_data_quality_analysis


	##RMA preprocessing
	#this have been done following oligo user guide
		#https://rdrr.io/bioc/oligo/f/inst/scripts/oug.pdf
	#Preprocessing refers to a series of complex statistical procedures applied to microarray data prior to dowstream analyses. These steps are required mainly for two reasons: A) technical artifacts are known to affect results, so background subtraction and normalization are used to minimize these issues; and B) there are MULTIPLE probes per probeset, therefore summarization to the probeset level is needed, so downstream analyses can be carried on
	
	#background subtraction
	#The oligo package implements background subtraction through the backgroundCorrect command. The method currently available is the one used in RMA, which treats the signal intensities of the perfect match (PM) probes as a convolution of noise and true signal. Additional methods will be available on future releases and choices will be made with the method argument (currently, the default is ’rma’).
		#see this link for PM probes (https://online.stat.psu.edu/stat555/node/28/)
	bg_AEset_human = backgroundCorrect(object=AEset_human, method="rma")
		#object: Object containing probe intensities to be preprocessed.
		#method: String determining which method to use at that preprocessing step

	#you can compare the probe intensities with and without background correction
	if(FALSE){ #in the oligo user guide, you can see how without background correction there is more noise
		par(mfcol(1,2))
		boxplot(bg_AEset_human, target="core")
		boxplot(AEset_human, target="core")
	}

	#Normalization
	#The normalize method provided by oligo allows the user to normalize the input data. Normalization is intended to remove from the intensity measures any systematic trends which arise from the microarray technology rather than from differences between the probes or between the target RNA samples hybridized to the arrays. Different normalization methods are available. The available options are given by normalizationMethods and the argument method in normalize is used to select the normalization approach to be used. 
	norm_bg_AEset_human <- normalize(object=bg_AEset_human)
		#object: A data object, typically containing microarray data.

	#Summarization
	#Summarizing the MULTIPLE probes per probeset to get a value at the probeset level. We use the same method used in RMA, which is median polish. 
	#summ_norm_bg_AEset_human = summarize(object=norm_bg_AEset_human, method="medianpolish") #NOT WORKING
		#object: Object containing probe intensities to be preprocessed.
		#method: String determining which method to use at that preprocessing step.

	#RMA vs MAS5
	#If you have HuGene arrays, you should use RMA (or GCRMA) normalization, since you don't have mismatch probes (although you have other types of control probes)
		#https://www.biostars.org/p/211389/

	#According the ArrayExpress vignette, we can use RMA normalisation to process the raw data from Affymetrix experiments.

	#In the literature, you will always be able to find examples where people state that one method performed better than another; here's an article extolling the virtues of MAS5. The important thing to remember is that they observed the improvement precisely once, under a specific set of conditions - you can't generalise to all cases from one good result. In general though, I disagree with your colleague: I'd say that RMA "is more often used nowadays. The essential differences between RMA and MAS5 are:
		#MAS5 normalises each array independently and sequentially; RMA as the name suggests (robust multi-array) uses a multi-chip model
		#MAS5 uses data from mismatch probes to calculate a "robust average", based on subtracting mismatch probe value from match probe value
		#RMA does not use the mismatch probes, because their intensities are often higher than the match probes, making them unreliable as indicators of non-specific binding
		#RMA values are in log2 units, MAS5 are not (so values are not directly comparable)
		#https://www.biostars.org/p/7687/

	#do all the steps of preprocessing with the RMA normalization using the oligo package
	AEsetnorm = oligo::rma(AEset_human, background=TRUE, normalize=TRUE, target="core")
		#the tutorial says affy:rma, but for me it does not work. I have to use oligo::rma. In addition, affy is much older.
			#https://support.bioconductor.org/p/71097/
		#oligo::rma
			#Robust Multichip Average preprocessing methodology. This strategy allows background subtraction, quantile normalization and summarization (via median-polish). In other words, during the RMA normalisation, a Tukey's 'median polish' is applied using information from all probesets.
			#object: Exon/HTA/Expression/Gene/SnpCnv-FeatureSet object.
			#background: Logical - perform RMA background correction?
			#normalize: Logical - perform quantile normalization?
			#target: Level of summarization (e.g., at the gene level... only for Exon/Gene arrays)
				#core would summarize to the gene level. However, if your array is an 'Exon' array, this may still only summarise to exon-level, requiring you to perform a further [manual] summarisation to gene-level.
					#we are going to do the manul summarization at the gene level anyways.
					#https://www.biostars.org/p/271379/#271588


	##filtering based on intensity
	#Once the RMA preprocessing is done, you can remove those probes with a very low expression. Microarray data commonly show a large number of probes in the background intensity range. These probes also do not change much across arrays. Hence they combine a low variance with a low intensity. Thus, they could end up being detected as differentially expressed although they are barely above the “detection” limit and are not very informative in general. A “soft” intensity based filtering is recommended by the limma (11,12) user guide in order to perform the differential expression analysis.
		#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#10_Filtering_based_on_intensity
	
	#we are going to create the histogram and then, we will filter if needed
	#calculate the median intensity per transcript across all arrays
	AEsetnorm_medians <- rowMedians(exprs(AEsetnorm))
	
	#plot the histogram
	pdf(paste(selected_ids_bat_studies, "/", selected_ids_bat_studies, "_hist_median_intensity.pdf", sep=""))
	hist_res <- hist(AEsetnorm_medians, 100, col = "cornsilk1", freq = FALSE, main = "Histogram of the median intensities", border = "antiquewhite4", xlab = "Median intensities")
	dev.off()

	#in case you need to filter after plotting the histogram, check this tutorial
		#https://www.bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#10_Filtering_based_on_intensity


	##filtering based on the experimental design
	#we are going to remove arrays that do not belong to BAT or have characteristics we are not interested in

	#create a vector with all the arrays before filtering
	selected_arrays = row.names(pData(AEsetnorm))

	#do the operations according to the dataset
	if(selected_ids_bat_studies == "E-GEOD-56635"){

		#select arrays that do not come from 
		bat_arrays_index = which(!grepl("White Adipocytes", AEsetnorm$Characteristics..cell.type., fixed=TRUE))
			#You can use pData to access to the phenotypic data (e.g., covariates) and meta-data (e.g., descriptions of covariates) associated with an experiment.
			#The quality metrics are better after removing WAT arrays.

		#create a phenotypic variable to define if the array comes from BAT or WAT
		#AEsetnorm$tissue = NA
		#AEsetnorm$tissue[bat_arrays_index] = "BAT"
		#AEsetnorm$tissue[-bat_arrays_index] = "WAT"

		#save the column used to filter
		#column_filter_names = "tissue"
		column_filter_names = "Characteristics..cell.type."

		#save the names of the selected arrays
		selected_arrays = selected_arrays[bat_arrays_index]
	}

	#select only the arrays we are interested
	AEsetnorm_filter = AEsetnorm[, selected_arrays]

	#check
	print("###############################################")
	print(paste(selected_ids_bat_studies, ": WE SELECTED THE INTEREST ARRAY?", sep=""))
	print(identical(colnames(exprs(AEsetnorm_filter)), selected_arrays))
	print("###############################################")


	##Quality assessment
	#To check the normalisation efficiency, we can run a quality assessment with arrayQualityMetrics. The arrayQualityMetrics package produces, through a single function call, a comprehensive HTML report of quality metrics about a microarray dataset [1, 2, 3]. The quality metrics are mainly on the per array level, i. e. they can be used to assess the relative quality of different arrays within a dataset. For example, different tissues.
	#It can be used for: i) assessing quality of a “raw” dataset, in order to get feedback on the experimental procedures that produced the data; (ii) assessing quality of a normalised dataset, in order to decide whether and how to use the dataset (or subsets of arrays in it) for subsequent data analysis.
	#Different types of microarray data (one colour, two colour, Affymetrix, Illumina) are represented by different object classes in Bioconductor. The function arrayQuality Metrics will work in the same way for all of them. Further information about its arguments can be found in its manual page.
		#https://bioconductor.org/packages/release/bioc/vignettes/arrayQualityMetrics/inst/doc/arrayQualityMetrics.pdf
	#Metadata about the arrays is shown at the top of the report as a table (see Figure 1). It is extracted from the phenoData slot of the data object supplied to arrayQuali tyMetrics. It can be useful to adjust the contents this slot before producing the report, and to make sure it contains the right quantity of information to make an informative report - not too much, not too little.
	
	#perform the quality assessment with and without the filter
	qanorm_filter = arrayQualityMetrics(expressionset=AEsetnorm_filter, outdir=paste(selected_ids_bat_studies, "/QAnorm_filter", sep=""), force=TRUE, do.logtransform = FALSE, intgroup=column_filter_names)
	qanorm_no_filter = arrayQualityMetrics(expressionset=AEsetnorm, outdir=paste(selected_ids_bat_studies, "/QAnorm_no_filter", sep=""), force=TRUE, do.logtransform = FALSE, intgroup=column_filter_names)
		#expressionset: a Bioconductor microarray dataset container. This can be an object of class ExpressionSet, AffyBatch, NChannelSet, ExpressionSetIllumina, RGList, MAList.
		#outdir: the name of the directory in which the report is created; a character of length 1.
		#force: if the directory named by ‘outdir’ already exists, then, if ‘force’ is ‘TRUE’, the directory is overwritten, otherwise an error is thrown; if the directory does not exist, the value of ‘force’ is irrelevant; a logical of length 1.
		#do.logtransform: indicates whether the data should be logarithm transformed before the analysis; a logical of length 1.
			#in our case, it has been already log transformed in the RMA normalization.
		#intgroup: the name of the sample covariate(s) used to draw a colour side bar next to the heatmap. The first element of ‘intgroup’ is also used define sample groups in other plots (boxplots, densities).  ‘intgroup’ should be a character vector, and its elements need to match the columns names of ‘phenoData(expressionset)’. If its length is 0, then the plots are not decorated with sample covariate information.
			#A useful feature of arrayQualityMetrics is the possibility to show the results in the context of an experimental factor of interest, i. e. a categorical variable associated with the arrays such as hybridisation date, treatment level or replicate number. Specifying a factor does not change how the quality metrics are computed. By setting the argument intgroup to contain the names of one or multiple columns of the data object’s phenoData slot5 , a bar on the side of the heatmap with colours representing the respective factors is added. Similarly, the colours of the boxplots and density plots reflect the levels of the first of the factors named by intgroup.

	#save the expression datasets with and without filter
	saveRDS(AEsetnorm_filter, paste(selected_ids_bat_studies, "/", selected_ids_bat_studies, "_expression_filter_1.Rds", sep=""))
	saveRDS(AEsetnorm, paste(selected_ids_bat_studies, "/", selected_ids_bat_studies, "_expression_non_filter.Rds", sep=""))
		#https://support.bioconductor.org/p/125920/
}
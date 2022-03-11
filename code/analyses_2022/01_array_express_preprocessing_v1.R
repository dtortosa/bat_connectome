#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



####################################################################################
########################## PREPROCESING GENE EXPRESSSION ###########################
####################################################################################

#Script for preprocesing gene expression data from ArrayExpress

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

#completely remove the previous folder with raw data
system(paste("rm -rf array_express_raw", sep=""))

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
sets[,which(colnames(sets) != "ExperimentFactors")]
	#we get the same, 15 experiments, just like we search in the web, the same IDs
	#https://www.ebi.ac.uk/arrayexpress/search.html?query=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22


##IMPORTANT: WE ARE SELECTING ONLY microarray studies, not RNA-seq studies
#There are much more studies based on microarray technology compare to RNA-seq. Therefore, we have more availability to select. In the case of RNA-seq, there are only 3 studies with raw data. In the case of microarray technology, there are much more studies. 
	#https://www.ebi.ac.uk/arrayexpress/search.html?query=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22

#Importantly, people do not recommend to combine different technologies, so we would need to do the analyses separately. For now, it makes sense to stick to microarray.
	#https://www.biostars.org/p/130986/

#If a reviewer asks for RNA-seq, you can use ArrayExpressHTS, which is a bioconductor package used for processing this type of data.
	#https://bioconductor.org/packages/release/bioc/html/ArrayExpressHTS.html


##select sets with raw data available
#Note that the results obtained from ArrayExpress do not consider RNA-seq datasets, so in all cases they are indicated as not having Raw nor Processed data.

#Among microarray, I prefer to select only those with raw data so I can process the data in the same way for all datasets. For example, if a dataset has been processed using RMA, expression has been log transformed, but not if the MAS5 processing was applied. Therefore gene expression of each study is not comparable. In our case, we are going to calculate tops in each dataset, but I think it is better to maintain the same treatment just in case...
sets_raw=sets[which(sets$Raw=="yes"),]
	#this considers only datasets with raw data STORED in ArrayExpress. Therefore, this discards cases with RNA-seq data stored in ENA. If a reviewer ask for it, we can try to include also this data.

#see
sets_raw[,which(colnames(sets_raw) != "ExperimentFactors")]

#check
print("#########################################")
print("CHECK THAT ALL DATASETS HAVE DATA FOR HUMANS AND INCLUDE RAW DATA")
print(length(which(sets_raw$Raw != "yes"))==0 & length(which(!grepl("Homo sapiens", sets_raw$Species)))==0)
print("#########################################")

#Note that ArrayExpress have stopped the regular imports of Gene Expression Omnibus (GEO) data into ArrayExpress. So, for GEOD datasets, the stored version may not be the latest version of the experiment.
	#we have to check that in each dataset


##select only those microarray experiments measuring gene expression in BAT-like tissues
#sets with BAT data
bat_datasets_info = list() #empty list to save info
#E-GEOD-56635
bat_datasets_info[[1]] = list(
	dataset="E-GEOD-56635",
	paper="Reprogrammed Functional Brown Adipocytes Ameliorate Insulin Resistance and Dyslipidemia in Diet-Induced Obesity and Type 2 Diabetes (doi:10.1016/j.stemcr.2015.08.007)",
	url="https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-56635/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=",
	details_bat="The original publication says that the accession number for the microarray data reported is GSE52817, but that number connect to data for the transformation of fibroblasts to osteoblasts. In contrast, the paper explains that created two types of brown adipocytes: i) iBA (induced brown adipocytes), coming from human induced pluripotent stem cells (iPSCs). After 12 days, 75% of the cells showed BA-like appearance (multilocular lipid droplets and abound mitochondria. They also showed high expression BAT markers (UCP1, CIDEA and DIO2) compared to human dermal fibroblasts. Stimulation of these cells with isoproterenol, a general agonist of β-adrenergic receptors, resulted in further elevation of expression level of UCP1 mRNA, which is a typical adrenergic response of BAs; dBA: Direct conversion from human fibroblasts. They used retroviruses to include reprogramming factor genes. CM-transduced genes were considered dBAs, being more than 90% of them UCP1 positive. They expressed mRNA of UCP1, CIDE and ADIPOQ at high levels compared to fibroblasts.  They also showed higher rate of oxygen consumption than fibroblasts. The oxygen uptake of dBAs was only partially suppressed by oligomycin, an inhibitor of mitochondrial ATP synthase, demonstrating a high rate of uncoupling respiration. Therefore, we can obtain gene expression data from cells showing BAT-like phenotypes: iBAs and dBAs. There are also mouse cells, but our script will select only human cells (see below). For human cells, RNA was obtained from WAs, iBAs, and dBAs, and after reverse-transcription, microarray analyses were performed using GeneChip human Gene 1.0 ST (Affymetrix) according to the manufacturer’s instruction.", 
	details_date="The dataset was updated 19 August 2015 for the last time in ArrayExpress, and the corresponding paper was published in 13 October 2015. Therefore, this is likely the last version of the dataset.")

#E-GEOD-19643
bat_datasets_info[[2]] = list(
	dataset="E-GEOD-19643",
	paper="PGC-1 alpha mediates differentiation of mesenchymal stem cells to brown adipose cells (doi:10.5551/jat.7401)",
	url="https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-19643/?query=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22",
	details_bat="They took mesenchymal stem cells from the bone marrow of a female donor (BMSCs). They induced adipocyte differentiation by treating BMSCs with a differentiation medium. They also performed osteogenic differentiation on the BMSCs. Before adipocyte differentiation, they transformed BMSCs with adenovirus containing green fluorescent protein (GFP) and human PGC-1 alpha cDNA or the antisense sequence of PGC-1 alpha, respectively. The adenovirus containing GFP (Ad-GFP) and adenovirus containing the antisense sequence of PGC-1 alpha (Ad-AS-PGC-1 alpha) were used as controls. Then, they used antibodies to detect several proteins: PGC-1 alpha, UCP1, UCP2, necdin, nuclear respiratory factor (NRF-1) and actin. They also measured mitochondrial mass, oxygen consumption and reactive oxygen species. Those BMSC cells infected with PGC-1 alpha viruses showed upregulation of expression of PGC-1 alpha mRNA. There is also an increase and mitochondrial mass and UCP1, along with a decrease in genes related to WAT. Therefore, these are BAT-like features. Note that this refers to bone marrow-derived MSCs infected with Ad-PGC-1α at a multiplicity of infection (m.o.i.) of 500 overnight. I understand this is the data we have (not 100 overnight which is also showed in the paper), because the information about this accession in NCBI indicates 500 overnight.",
	details_date="The dataset was updated 10 June 2011 for the last time in ArrayExpress, and the corresponding paper was published in 5 August 2011. Therefore, this is likely the last version of the dataset.")

#E-GEOD-27657
bat_datasets_info[[3]] = list(dataset="E-GEOD-27657",
	paper="Gene expression in human brown adipose tissue (doi:10.3892/ijmm.2010.566)",
	url="https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-27657/",
	details_bat="Biopses of adipose tissue of brown-red color were obtained from an area close to the isthmus region of the thyroid gland. In each patient, biopses were also taken from the subcutaneous depot in the surgical incision area and were considered as WAT. They used staining and also antibodies against UCP1 or CKMT1. Gene expression was evaluated using the Human Genome U133Plus2.0 DNA microarray (Affymetrics). UCP1 was undetectable in all WAT samples, but it was detectable in half of the brown-red biopsies. Biopsies were classified as BAT if UCP1 expression was >0.5 relative to the reference gene RPLP0. They analyzed genes deferentially expressed in BAT and WAT. They found that UCP1 gene was overexpressed in sample classified as BAT. Other 16 genes were also expressed at higher levels in BAT compared to WAT. According to Gene Ontology (GO) analysis, five of these genes encode proteins that have their main function in mitochondria. In the abstract, they say 'Microarray analysis of 9 paired human BAT and WAT samples showed that 17 genes had at least a 4-fold higher expression in BAT compared to WAT and five of them'. Therefore, I understand that these are the 18 samples I have in the expression file, and those belonging to perthyroid are those biopsies showing BAT features.",
	details_date="The dataset was updated 2 June 2014 for the last time in ArrayExpress, and the two related papers were published in 1 December 2010 and 1 February 2013. Therefore, this is likely the last version of the dataset.")

#E-GEOD-54280
bat_datasets_info[[4]] = list(dataset="E-GEOD-54280", 
	paper="Comparative gene array analysis of progenitor cells from human paired deep neck and subcutaneous adipose tissue (doi: 10.1016/j.mce.2014.07.011)",
	url="https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-54280/",
	details_bat="They took deep neck adipose tissue from the region surrounding the carotid sheath. From the same patient, subcutaneous neck adipose tissue was obtained from the neck at the surgical incision. BAT was defined by the presence of multivacuolar adipocytes staining positive for UCP1. Progenitor cells were isolated by collagenase digestion and subsequently transferred to cell culture. After reaching subconfluency (7 days), they harvested RNA and analyzed differences in gene expression by microarray analysis. The total RNAs were amplified and labeled following the Whole Transcript (WT) Sense Target Labeling Assay (http://www.affymetrix.com). Labeled ssDNA was hybridized to Human Gene 1.0 ST Affymetrix GeneChip arrays (Affymetrix). The regions of deep neck adipose tissue that contained multivacuolar adipocytes stained positively for UCP1, but no staining was observed in univacuolar adipocytes in deep neck and subcutaneous neck tissue. They found an enrichment of UCP1 and PRDM16 mRNA expression in deep neck tissue samples compared with subcutaneous, but high variance. The same was found for LHX8, a BAT marker. However, they found a higher significant expression of UCP1 and PRDM16 in adipocytes derived from deep neck cells after differentiation (they were exposed to a differentiation medium) compared with adipocytes derved from subcutaneous cells. Same for RDM16 and LHX8. Leptin, which is a WAT marker, was more expressed in adipocytes derived from subcutaneous cells compared to adipocytes derived from the neck. Although there were not great functional differences respect oxygen consumption and proton leak. Finally, they compared gene expression in microarray between progenitor cells of deep and subcutaneous neck in six patients. These are our patients, I understand that they took only 6 of the 12 individuals and then took 2 samples for each one (deep vs. subcutaneous neck), giving our 12 samples. Genes previously associated with BAT or browning, like BMP4, EBF2, FABP3 and PDGFR alpha were higher expressed in progenitor cells from the deep neck adipose depot compared with cells from the subcutaneous neck tissue. In the same line, genes which have been described as marker gene for white adipose tissue as well as genes implicated in the determination of white adipocytes including VDR, HOXC8, SHOX2, TWIST1 and PAX3 were higher expressed in subcutaneous neck compared with deep neck progenitor cells. Therefore, deep neck samples seems to be closer to BAT compared to subcutaneous neck samples.",
	details_date="The dataset was updated 8 August 2014 for the last time in ArrayExpress, and the corresponding paper was published in September 2014. Therefore, this is likely the last version of the dataset.")

#create a vector with the IDs
ids_bat_studies = as.vector(unlist(sapply(X=bat_datasets_info, "[", "dataset")))

#save the list with information
saveRDS(object=bat_datasets_info, file="info_datasets.Rds")

#IMPORTANT: Files of some studies cannot be loaded, the error is duplicated row.names. Maybe you can dowload separately each part (getAE) and create the ExpressionSet file by yourself
	#E-MTAB-25 fail, Error in `.rowNamesDF<-`(x, value = value) :  duplicate 'row.names' are not allowed
		#maybe: https://www.biostars.org/p/62988/



#######################################
#### OBTAIN DATA FOR EACH DATASET #####
#######################################

##write a function to do that
#selected_ids_bat_studies=ids_bat_studies[1] #for debugging
extract_expression_data = function(selected_ids_bat_studies){

	#extend the time window during which a function can download, for example, download.file, which is used by ArrayExpress function. We need this inside the function, so the same timeout is applied for each independent run.
	getOption('timeout')
	options(timeout=1000)
	getOption('timeout')
		#https://stackoverflow.com/questions/35282928/how-do-i-set-a-timeout-for-utilsdownload-file-in-r

	#make a directory for the selected dataset
	system(paste("mkdir -p ", selected_ids_bat_studies, sep=""))
		#use mkdir -p because if you don't and run again, you will get an error that it exists. "p" flags is for "no error if existing, make parent directories as needed"

	#see one specific dataset
	print("###############################################")
	print(paste("INFO", selected_ids_bat_studies, ": ", sep="")); print(sets_raw[which(sets_raw$ID == selected_ids_bat_studies), 1:7])
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
	if(!grepl("hugene", AEset_human@annotation, fixed=TRUE) & !grepl("hg", AEset_human@annotation, fixed=TRUE)){
		print("###############################################")
		stop(paste(selected_ids_bat_studies, " HAS NO HUMAN GENES", sep=""))
		print("###############################################")
	} #According to some people, if you have HuGene arrays, you should use RMA (or GCRMA) normalization, since you don't have mismatch probes (although you have other types of control probes).
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
	if(FALSE){
		bg_AEset_human = backgroundCorrect(object=AEset_human, method="rma")
			#object: Object containing probe intensities to be preprocessed.
			#method: String determining which method to use at that preprocessing step
	}

	#you can compare the probe intensities with and without background correction
	if(FALSE){ #in the oligo user guide, you can see how without background correction there is more noise
		par(mfcol(1,2))
		boxplot(bg_AEset_human, target="core")
		boxplot(AEset_human, target="core")
	}

	#Normalization
	#The normalize method provided by oligo allows the user to normalize the input data. Normalization is intended to remove from the intensity measures any systematic trends which arise from the microarray technology rather than from differences between the probes or between the target RNA samples hybridized to the arrays. Different normalization methods are available. The available options are given by normalizationMethods and the argument method in normalize is used to select the normalization approach to be used.
	if(FALSE){
		norm_bg_AEset_human <- normalize(object=bg_AEset_human)
			#object: A data object, typically containing microarray data.
	}
	
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
	AEsetnorm = oligo::rma(AEset_human, background=TRUE, normalize=TRUE)
		#the tutorial says affy:rma, but for me it does not work. I have to use oligo::rma. In addition, affy is much older.
			#https://support.bioconductor.org/p/71097/
		#oligo::rma
			#Robust Multichip Average preprocessing methodology. This strategy allows background subtraction, quantile normalization and summarization (via median-polish). In other words, during the RMA normalisation, a Tukey's 'median polish' is applied using information from all probesets.
			#object: Exon/HTA/Expression/Gene/SnpCnv-FeatureSet object.
			#background: Logical - perform RMA background correction?
			#normalize: Logical - perform quantile normalization?
			#target: Level of summarization (e.g., at the gene level... only for Exon/Gene arrays)
				#core would summarize to the gene level. However, if your array is an 'Exon' array, this may still only summarise to exon-level, requiring you to perform a further [manual] summarisation to gene-level.
					#the default is "core"
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

		#select arrays that do not come from white adypocytes 
		bat_arrays_index = which(!grepl("White Adipocytes", AEsetnorm$Characteristics..cell.type., fixed=TRUE))
			#You can use pData to access to the phenotypic data (e.g., covariates) and meta-data (e.g., descriptions of covariates) associated with an experiment.

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
	if(selected_ids_bat_studies == "E-GEOD-19643"){

		#select arrays for PGC-1 alpha-transformed cells
		bat_arrays_index = which(grepl("PGC-1 alpha", AEsetnorm$Comment..Sample_source_name., fixed=TRUE))
			#PGC-1 alpha transformed mesenchymal cells are those showing BAT-like features.
			#You can use pData to access to the phenotypic data (e.g., covariates) and meta-data (e.g., descriptions of covariates) associated with an experiment.

		#notes after quality check for this dataset
		#in general, there are not great signals of problem in specific arrays. The processing worked well.
		#note, however, that the two arrays of PGC-1 alpha, which shows BAT-like features, do not constitute a clear group opposite to the non-PGC-1 alpha arrays, i.e., controls. This is something to consider, because maybe these cells are not different enough from the controls. Therefore, this is another reason to the the analyses considering only BAT biopsies studies.

		#save the column used to filter
		column_filter_names = "Comment..Sample_source_name."

		#save the names of the selected arrays
		selected_arrays = selected_arrays[bat_arrays_index]
	}
	if(selected_ids_bat_studies == "E-GEOD-27657"){

		#select arrays for BAT
		bat_arrays_index = which(grepl("perithyroid", AEsetnorm$Characteristics..organism.part., fixed=TRUE))
			#This includes female and male subjects
			#You can use pData to access to the phenotypic data (e.g., covariates) and meta-data (e.g., descriptions of covariates) associated with an experiment.

		#save the column used to filter
		column_filter_names = "Characteristics..organism.part."

		#save the names of the selected arrays
		selected_arrays = selected_arrays[bat_arrays_index]

		#remove three individuals with perithyroid data which seems to be an outlier
		selected_arrays = selected_arrays[which(!selected_arrays %in% c("GSM685080.CEL", "GSM685075.CEL", "GSM685074.CEL"))]
			#when analyzing all arrays (including WAT), these arrays are considered as outliers
				#GSM685080.CEL according to the distances between arrays.
				#GSM685075.CEL according to the distribution of the intensity.
				#GSM685074.CEL according to the distances between arrays and the distribution of the intensity.
			#the removal of these arrays does not qualitatively change the results.
	}
	if(selected_ids_bat_studies == "E-GEOD-54280"){

		#select arrays for BAT
		bat_arrays_index = which(grepl("deep neck", AEsetnorm$Comment..Sample_source_name., fixed=TRUE))
			#This includes male and female subjects
			#You can use pData to access to the phenotypic data (e.g., covariates) and meta-data (e.g., descriptions of covariates) associated with an experiment.

		#save the column used to filter
		column_filter_names = "Comment..Sample_source_name."

		#save the names of the selected arrays
		selected_arrays = selected_arrays[bat_arrays_index]

		#remove a sample that it is very close to WAT samples
		selected_arrays = selected_arrays[which(!selected_arrays %in% c("GSM1311783_01-11_P1.CEL"))]
			#GSM1311783_01-11_P1.CEL is close to WAT arrays according to the PCA. After its removal, BAT and WAT sample are more clearly differentiated. Note that GSM1311787_23-12_P1.CEL is considered now an outlier according to the distribution of intensity, but this is only when BAT sample are considered in isolation, when WAT is considered, there are no outliers. In addition, if you see the distribution of intensity in density plots with only BAT, there is not array that differentiate a lot from the rest, so we are going to end here.
	}

	#select only the arrays we are interested
	AEsetnorm_filter = AEsetnorm[, selected_arrays]

	#check
	print("###############################################")
	print(paste(selected_ids_bat_studies, ": WE SELECTED THE INTEREST ARRAY?", sep="")); print(identical(colnames(exprs(AEsetnorm_filter)), selected_arrays))
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

	#return the name of the dataset processed
	return(selected_ids_bat_studies)
}



#####################################
###### PARALLELIZE THE PROCESS ######
#####################################

#BAT studies
ids_bat_studies

#set up cluster
clust <- makeCluster(length(ids_bat_studies), outfile="") #outfile let you to see the output in the terminal "https://blog.revolutionanalytics.com/2015/02/monitoring-progress-of-a-foreach-parallel-job.html"
registerDoParallel(clust)

#run the function for all populations
ids_bat_studies_processed=foreach(i=ids_bat_studies, .packages=c("ArrayExpress", "arrayQualityMetrics")) %dopar% {
    extract_expression_data(selected_ids_bat_studies=i)
}

#stop the cluster 
stopCluster(clust)

#check
print("###############################################")
print(paste("WE PROCESSED ALL SELECTED DATASETS?", sep="")); print(summary(ids_bat_studies == ids_bat_studies_processed))
print("###############################################")
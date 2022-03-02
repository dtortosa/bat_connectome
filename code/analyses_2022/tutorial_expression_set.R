#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



#######################################################################
########################## BIOBASE TUTORIAL ###########################
#######################################################################

#tutorial about using ExpressionSet objects
	#https://www.bioconductor.org/packages/devel/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf



###########################
#### REQUIRES PACKAGES ####
###########################

require(Biobase)



##############################################
#### BUILDING EXPRESSION SET FROM SCRATCH ####
##############################################

#### Assay Data ####
#A likely scenario is that your assay data is in a tab-delimited text file (as exported from a spreadsheet, for instance) with rows corresponding to features and columns to samples. The strategy is to read this file into R using the read.table command, converting the result to a matrix
dataDirectory <- system.file("extdata", package="Biobase") #get the path where the example data of biobase is located
exprsFile <- file.path(dataDirectory, "exprsData.txt") #get the complete path with the file name
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE)) #load the file

#take a look
str(exprs)
class(exprs)
dim(exprs)
colnames(exprs)
head(exprs[,1:5])
	#I guess the rows are the genes and the columns are the individuals, because some rows have as name IL...

#At this point, we can create a minimal ExpressionSet object using the ExpressionSet constructor:
minimalSet <- ExpressionSet(assayData=exprs)

#We’ll get more benefit from expression sets by creating a richer object that coordinates phenotypic and other data with our expression data, as outlined in the following sections.



#### Phenotypic data ####

#Phenotypic data summarizes information about the samples (e.g., sex, age, and treatment status; referred to as ‘covariates’). The information describing the samples can be represented as a table with S rows and V columns, where V is the number of covariates.
pDataFile <- file.path(dataDirectory, "pData.txt")
pData <- read.table(pDataFile, row.names=1, header=TRUE, sep="\t")
dim(pData)
rownames(pData)
colnames(pData)
summary(pData)

#There are three columns of data, and 26 rows. Note that the number of rows of phenotypic data match the number of columns of expression data, and indeed that the row and column names are identically ordered:
all(rownames(pData)==colnames(exprs))

#This is an essential feature of the relationship between the assay and phenotype data; ExpressionSet will complain if these names do not match.

#Phenotypic data can take on a number of different forms. For instance, some covariates might reasonably be represented as numeric values. Other covariates (e.g., gender, tissue type, or cancer status) might better be represented as factor objects (see the help page for factor for more information). It is especially important that the phenotypic data are encoded correctly; the colClasses argument to read.table can be helpful in correctly inputing (and ignoring, if desired) columns from the file.
sapply(pData, class)
pData[c(15,20), c("gender", "type")]
pData[pData$score>0.8,]

#Investigators often find that the meaning of simple column names does not provide enough information about the covariate – What is the cryptic name supposed to represent? What units are the covariates measured in? We can create a data frame containing such meta-data (or read the information from a file using read.table)
metadata <- data.frame(labelDescription= c("Patient gender", "Case/control status", "Tumor progress on XYZ scale"), row.names=c("gender", "type", "score"))
	#This creates a data.frame object with a single column called labelDescription, and with row names identical to the column names of the data.frame containing the phenotypic data. The column labelDescription must be present; other columns are optional.

#Bioconductor’s Biobase package provides a class called AnnotatedDataFrame that conveniently stores and manipulates the phenotypic data and its metadata in a coordinated fashion
phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
phenoData

#Some useful operations on an AnnotatedDataFrame include sampleNames, pData (to extract the original pData data.frame), and varMetadata. In addition, AnnotatedDataFrame objects can be subset much like a data.frame:
head(pData(phenoData))
phenoData[c("A","Z"),"gender"]
pData(phenoData[phenoData$score>0.8,])


#### Annotations and feature data ####

#Meta-data on features is as important as meta-data on samples, and can be very large and diverse. A single chip design (i.e., collection of features) is likely to be used in many different experiments, and it would be inefficient to repeatedly collect and coordinate the same metadata for each ExpressionSet instance. Instead, the ideas is to construct specialized meta-data packages for each type of chip or instrument. Many of these packages are available from the Bioconductor web site. These packages contain information such as the gene name, symbol and chromosomal location. There are other meta-data packages that contain the information that is provided by other initiatives such as GO and KEGG. The annotate and AnnotationDbi packages provides basic data manipulation tools for the meta-data packages.

#The appropriate way to create annotation data for features is very straight-forward: we provide a character string identifying the type of chip used in the experiment. For instance, the data we are using is from the Affymetrix hgu95av2 chip
annotation <- "hgu95av2"

#It is also possible to record information about features that are unique to the experiment (e.g., flagging particularly relevant features). This is done by creating or modifying an AnnotatedDataFrame like that for phenoData but with row names of the AnnotatedDataFrame matching rows of the assay data.


#### Experiment description ####

#Basic description about the experiment (e.g., the investigator or lab where the experiment was done, an overall title, and other notes) can be recorded by creating a MIAME object. One way to create a MIAME object is to use the new function:
experimentData <- new("MIAME", name="Pierre Fermat", lab="Francis Galton Lab", contact="pfermat@lab.not.exist", title="Smoking-Cancer Experiment", abstract="An example ExpressionSet", url="www.lab.not.exist", other=list( notes="Created from text files" ))
	#Usually, new takes as arguments the class name and pairs of names and values corresponding to different slots in the class; consult the help page for MIAME for details of available slots.



#### Assembling an ExpressionSet ####

#An ExpressionSet object is created by assembling its component parts and callng the ExpressionSet constructor:
exampleSet <- ExpressionSet(assayData=exprs, phenoData=phenoData, experimentData=experimentData, annotation=annotation)
	#Note that the names on the right of each equal sign can refer to any object of appropriate class for the argument. See the help page for ExpressionSet for more information.

#We created a rich data object to coordinate diverse sources of information. Less rich objects can be created by providing less information. As mentioned earlier, a minimal expression set can be created with
minimalSet <- ExpressionSet(assayData=exprs)

#Of course this object has no information about phenotypic or feature data, or about the chip used for the assay.



###############################
#### EXPRESSION SET BASICS ####
###############################

#Now that you have an ExpressionSet instance, let’s explore some of the basic operations. You can get an overview of the structure and available methods for ExpressionSet objects by reading the help page
help("ExpressionSet-class")

#When you print an ExpressionSet object, a brief summary of the contents of the object is displayed (displaying the entire object would fill your screen with numbers):
exampleSet


#### Accessing Data Elements ####
#A number of accessor functions are available to extract data from an ExpressionSet instance. You can access the columns of the phenotype data (an AnnotatedDataFrame instance) using $
exampleSet$gender[1:5]
exampleSet$gender[1:5] == "Female"

#You can retrieve the names of the features using featureNames. For many microarray datasets, the feature names are the probe set identifiers.
featureNames(exampleSet)[1:5] #genes

#The unique identifiers of the samples in the data set are available via the sampleNames method.
sampleNames(exampleSet)[1:5] #study subjects

#The varLabels method lists the column names of the phenotype data:
varLabels(exampleSet)

#Extract the expression matrix of sample information using exprs:
mat <- exprs(exampleSet)
dim(mat)

#Probably the most useful operation to perform on ExpressionSet objects is subsetting. Subsetting an ExpressionSet is very similar to subsetting the expression matrix that is contained within the ExpressionSet, the first argument subsets the features and the second argument subsets the samples. Here are some examples: Create a new ExpressionSet consisting of the 5 features and the first 3 samples
vv <- exampleSet[1:5, 1:3]
dim(vv)
featureNames(vv)
sampleNames(vv)

#Create a subset consisting of only the male samples:
males <- exampleSet[, exampleSet$gender == "Male"]
males
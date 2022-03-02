#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



##############################################################################
########################## ANNOTATION BIOCONDUCTOR ###########################
##############################################################################

#tutorial about using ExpressionSet objects
	#https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf

#The most popular annotation packages have been modified so that they can make use of a new set of methods to more easily access their contents. These four methods are named: columns, keytypes, keys and select. And they are described in this vignette. They can currently be used with all chip, organism, and TxDb packages along with the popular GO.db package

#For the older less popular packages, there are still conventient ways to retrieve the data. The How to use bimaps from the ".db" annotation packages vignette in the AnnotationDbi package is a key reference for learnign about how to use bimap objects.

#Finally, all of the ‘.db’ (and most other Bioconductor annotation packages) are updated every 6 months corresponding to each release of Bioconductor. Exceptions are made for packages where the actual resources that the packages are based on have not themselves been updated



###########################
#### REQUIRES PACKAGES ####
###########################

require(hgu95av2.db)



####################################################
#### ANNOTATIONDB OBJECTS AND THE SELECT METHOD ####
####################################################

#As previously mentioned, a new set of methods have been added that allow a simpler way of extracting identifier based annotations. All the annotation packages that support these new methods expose an object named exactly the same as the package itself. These objects are collectively called AnntoationDb objects for the class that they all inherit from. The more specific classes (the ones that you will actually see in the wild) have names like OrgDb, ChipDb or TxDb objects. These names correspond to the kind of package (and underlying schema) being represented. The methods that can be applied to all of these objects are columns, keys, keytypes and select.

#In addition, another accessor has recently been added which allows extraction of one column at at time. the mapIds method allows users to extract data into either a named character vector, a list or even a SimpleCharacterList. This method should work with all the different kinds of AnntoationDb objects described below.



##############################################
#### CHIPDB OBJECTS AND THE SELECT METHOD ####
##############################################

#An extremely common kind of Annotation package is the so called platform based or chip based package type. This package is intended to make the manufacturer labels for a series of probes or probesets to a wide range of gene-based features. A package of this kind will load an ChipDb object. Below is a set of examples to show how you might use the standard 4 methods to interact with an object of this type.
	#this is our case

#First we need to load the package.
require(hgu95av2.db) 

#If we list the contents of this package, we can see that one of the many things loaded is an object named after the package "hgu95av2.db":
ls("package:hgu95av2.db")

#We can look at this object to learn more about it:
hgu95av2.db

#If we want to know what kinds of data are retriveable via select, then we should use the columns method like this
columns(hgu95av2.db)

#If we are further curious to know more about those values for columns, we can consult the help pages. Asking about any of these values will pull up a manual page describing the different fields and what they mean.
help("SYMBOL")

#If we are curious about what kinds of fields we could potentiall use as keys to query the database, we can use the keytypes method. In a perfect world, this method will return values very similar to what was returned by columns, but in reality, some kinds of values make poor keys and so this list is often shorter
keytypes(hgu95av2.db)

#If we want to extract some sample keys of a particular type, we can use the keys method.
head(keys(hgu95av2.db, keytype="SYMBOL"))

#And finally, if we have some keys, we can use select to extract them. By simply using appropriate argument values with select we can specify what keys we want to look up values for (keys), what we want returned back (columns) and the type of keys that we are passing in (keytype)
#1st get some example keys
k <- head(keys(hgu95av2.db, keytype="PROBEID"))
# then call select
select(hgu95av2.db, keys=k, columns=c("SYMBOL","GENENAME"), keytype="PROBEID")
	#And as you can see, when you call the code above, select will try to return a data.frame with all the things you asked for matched up to each other. 

#Finally if you wanted to extract only one column of data you could instead use the mapIds method like this:
#1st get some example keys
k <- head(keys(hgu95av2.db, keytype="PROBEID"))
# then call mapIds
mapIds(hgu95av2.db, keys=k, column=c("GENENAME"), keytype="PROBEID")



#############################################
#### ORGDB OBJECTS AND THE SELECT METHOD ####
#############################################

#An organism level package (an ‘org’ package) uses a central gene identifier (e.g. Entrez Gene id) and contains mappings between this identifier and other kinds of identifiers (e.g. GenBank or Uniprot accession number, RefSeq id, etc.). The name of an org package is always of the form org.<Ab>.<id>.db (e.g. org.Sc.sgd.db) where <Ab> is a 2-letter abbreviation of the organism (e.g. Sc for Saccharomyces cerevisiae) and <id> is an abbreviation (in lower-case) describing the type of central identifier (e.g. sgd for gene identifiers assigned by the Saccharomyces Genome Database, or eg for Entrez Gene ids).

#Just as the chip packages load a ChipDb object, the org packages will load a OrgDb object. The following exercise should acquaint you with the use of these methods in the context of an organism package.

#Display the OrgDb object for the org.Hs.eg.db package. 
library(org.Hs.eg.db)

#Use the columns method to discover which sorts of annotations can be extracted from it. Is this the same as the result from the keytypes method? Use the keytypes method to find out.
columns(org.Hs.eg.db)
help("SYMBOL") #for explanation of these columns and keytypes values

#Finally, use the keys method to extract UNIPROT identifiers and then pass those keys in to the select method in such a way that you extract the gene symbol and KEGG pathway information for each. Use the help system as needed to learn which values to pass in to columns in order to achieve this.
keytypes(org.Hs.eg.db)
uniKeys <- head(keys(org.Hs.eg.db, keytype="UNIPROT"))
cols <- c("SYMBOL", "PATH")
	#PATH: KEGG Pathway Identifiers
	#SYMBOL: The official gene symbol
select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="UNIPROT")

#So how could you use select to annotate your results? This next exercise should help you to understand how that should generally work.
#Please run the following code snippet (which will load a fake data result that I have provided for the purposes of illustration):
load(system.file("extdata", "resultTable.Rda", package="AnnotationDbi"))
head(resultTable)

#The rownames of this table happen to provide entrez gene identifiers for each row (for human), having each one a p-value just like a real dataset with significant differences in expression between treatments, tissues... Find the gene symbol and gene name for each of the rows in resultTable
annots <- select(org.Hs.eg.db, keys=rownames(resultTable), columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
#then use the merge method to attach those annotations to it
resultTable <- merge(resultTable, annots, by.x=0, by.y="ENTREZID")
head(resultTable)



#################################
#### USING SELECT WITH GO.DB ####
#################################

#When you load the GO.db package, a GODb object is also loaded. This allows you to use the columns, keys, keytypes and select methods on the contents of the GO ontology. So if for example, you had a few GO IDs and wanted to know more about it, you could do it like... STOPPED HERE.



#########################################
#### USING SELECT WITH TXDB PACKAGES ####
#########################################



##########################################
#### USING SELECT WITH ENSDB PACKAGES ####
##########################################
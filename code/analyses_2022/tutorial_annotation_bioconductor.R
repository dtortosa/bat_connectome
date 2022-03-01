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

require(pd.hugene.1.0.st.v1)



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

#First we need to load the package. We are going to use the annotation package of one of the datasets in ArrayExpress we are going to use for the BAT connectome (pd.hugene.1.0.st.v1). This the annotation of E-GEOD-54280 dataset.
	#http://www.bioconductor.org/packages/release/data/annotation/html/pd.hugene.1.0.st.v1.html
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-67297/?query=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22

#If we list the contents of this package, we can see that one of the many things loaded is an object named after the package "hgu95av2.db":
ls("package:pd.hugene.1.0.st.v1")

#ERROR
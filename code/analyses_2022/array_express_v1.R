
#tutorial: https://www.biostars.org/p/123233/
#Read this paper to get the sense of how to import an ArrayExpress dataset into R. This manual to know what an NChannel object is. Also, read Section 8 of this manual on how to identify differentially expressed genes from an ArrayExpress dataset.

require(ArrayExpress)
require(arrayQualityMetrics)
require(biomaRt)

#tutorial ArrayExpress
	#ArrayExpress general tutorial about the server
		#https://www.ebi.ac.uk/training/online/courses/array-express-discover-functional-genomics-data-quickly-and-easily/#vf-tabs__section--contents
	#ArrayExpress package
		#http://www.bioconductor.org/packages/release/bioc/vignettes/ArrayExpress/inst/doc/ArrayExpress.pdf

#set working directory to results, because ArrayExpress save stuff during the analyses automatically
setwd("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/results_2022")

#Use ArrayExpress for loading the data:
#If the identifier refers to an Affymetrix experiment, the output is an AffyBatch, if it refers to a one-colour experiment using a platform other than Affymetrix, the output is an ExpressionSet. The ArrayExpress function extracts feature intensity summaries from columns of the raw data files based on the common conventions for the data file sources. If the data source is not recognized, or the file does not have the expected column names, the user is asked to explicitly provide the name of the column(s) to extract, for instance, ‘Cy3 Median’. In some cases, there is a mismatch between the sample or feature annotations and the intensity data files; in such cases, a warning is emitted, the phenoData and/or featureData components are left empty and an incomplete (but syntactically valid) object is returned. 


#make a query for BAT
sets=queryAE(keywords="brown+adipose+tissue", species="homo+sapiens")

#take a look
str(sets)
	#we get the same, 15 experiments, just like we search in the web, the same IDs
	#https://www.ebi.ac.uk/arrayexpress/search.html?query=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22

#see one specific dataset
sets[which(sets$ID == "E-GEOD-54280"),]
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-67297/?query=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22

#extend the time window during which a function can download, for example, download.file, which is used by ArrayExpress function
getOption('timeout')
options(timeout=500)
getOption('timeout')
	#https://stackoverflow.com/questions/35282928/how-do-i-set-a-timeout-for-utilsdownload-file-in-r

#get this dataset
AEset=ArrayExpress("E-GEOD-54280")

#we have a list of two elements? two species?
#names(AEset) 
	#"A-AFFY-130" and "A-AFFY-141". Dataset for mouse and human respectively
		#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-54280/samples/?keywords=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=
#AEset_human = AEset$"A-AFFY-141"
AEset_human = AEset

#see manufacturer
AEset_human@manufacturer

#A-AFFY-141 is an Affymetrix experiment, so according the tutorial, we can use RMA normalisation to process the raw data, please read the rma function help.
#require(affy)
AEsetnorm = oligo::rma(AEset_human)
	#the tutorial says affy:rma, but for me it does not work. I have to use oligo::rma.
		#https://support.bioconductor.org/p/71097/

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


expression_matrix = exprs(AEsetnorm)
head(expression_matrix)

#we have a value of expression per subject, so we have to find a way to summarize and get a value across all subject per gene.
	#it makes sense to average expression across sex? but what about case/controls? or cold vs non cold?

average_gene = data.frame(apply(X=expression_matrix, MARGIN=1, FUN=mean))

#we have to convert the probs ID of affy_hugene_1_0_st_v1 to the gene symbols used by Yuval
	#https://www.biostars.org/p/69597/
	#https://support.bioconductor.org/p/42839/

require(biomaRt)
# replace the affyID with gene symbol
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",host = "www.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")

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


##DO THIS TUTORIAL ABOUT ANNOTATION!
	#https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
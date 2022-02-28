
#tutorial: https://www.biostars.org/p/123233/
#Read this paper to get the sense of how to import an ArrayExpress dataset into R. This manual to know what an NChannel object is. Also, read Section 8 of this manual on how to identify differentially expressed genes from an ArrayExpress dataset.
	#https://watermark.silverchair.com/btp354.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAuIwggLeBgkqhkiG9w0BBwagggLPMIICywIBADCCAsQGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQM_mrmYDgbljRfijwhAgEQgIIClYKZOVVDTEf0XXhbwdTLQnK9-7-3bhpFsi3C0P7TTC8RRjzO7nuXQr3Fw2r0ekyK9cwyi2uguv7FSQU2vKGM34UGpO65CnYjiibNzNnSD9UAyvAYyu1LBifna_uu15sWalDlsW-Vd_ILCXiF3GYE52YcUB4Es35STYUJ-QFGWdK34cQ5sFmfXLIeeODljmptq998MhBJ9yQE4QD_I0s3Y9gSNRz2Cq2pprXS3U635FYzH5SDXQPnKyIjzkR4-Ga0EuiTR3c0a4BZmaj3y56Qz7lLPoWkqpSOOBjUKsqS4REH-s5jnZyBLi5NYYxq1VefI1N6wdF-Fc2p8_x_XBkA3PO9PCSeJDZ1IwZpXiW0qYIu76HS0tJAFIAXKELdUynqpAmuBHjF3pmIOZzpCbYzGvQWE7sKlvYwe_RFoRIyCnitiKnK4gQQ3xgo7QXZTZJQavGsLfuEp2zIWlrKr3NTgWA388nmeGS-LtC3YfPTw6Ex0M_hKQR-pTzv2Vyz9jA7ajzbuGoLiuM--wmrWUJg21Za_by8r1xzKOt8M8woEekr2OlQLuaVA6ZmejUwSx3dLSF80FV6MQ0WiOI8SWI_02xEn0hCh_o9Ltdm_bYh9_4T_5RM3GTOmVSe_CEIJkOh9O8cDbk0kgQjBLZ61FWk1-O-cwWZX5dNpjcQ28GVFVCHtmPoUExwmez7_4-Jq-m2_Y-FyJb60TgDacqHTigZ9UFUdonUyPyHDQL5CdCv_R3ujMF6csdPht01NLJPTxGVUSKBhM_bfjl-gvzZ20CVTPWjtCpV5Ao3Oqpcb4jaJtMR7GXa6kPdq10CrU7thsojpjy7ljZPdlqgAeJRsf0wBRnr_9N713ACuB2Sp3eY8rGry6WvKEE
	#http://www.bioconductor.org/packages/release/bioc/manuals/Biobase/man/Biobase.pdf
	#http://www.bioconductor.org/packages/release/bioc/vignettes/ArrayExpress/inst/doc/ArrayExpress.pdf

#Use ArrayExpress for loading the data:
#As E-ATMX-18 is a two-colour experiment, the returned R object is of class NChannelSet. If the identifier refers to an Affymetrix experiment, the output is an AffyBatch, if it refers to a one-colour experiment using a platform other than Affymetrix, the output is an ExpressionSet. The ArrayExpress function extracts feature intensity summaries from columns of the raw data files based on the common conventions for the data file sources. If the data source is not recognized, or the file does not have the expected column names, the user is asked to explicitly provide the name of the column(s) to extract, for instance, ‘Cy3 Median’. In some cases, there is a mismatch between the sample or feature annotations and the intensity data files; in such cases, a warning is emitted, the phenoData and/or featureData components are left empty and an incomplete (but syntactically valid) object is returned. 
AEset = ArrayExpress("E-ATMX-18")


#make a query for BAT
sets=queryAE(keywords="brown+adipose+tissue", species="homo+sapiens")

#take a look
str(sets)
	#we get the same, 15 experiments, just like we search in the web, the same IDs
	#https://www.ebi.ac.uk/arrayexpress/search.html?query=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22

#see one specific dataset
sets[which(sets$ID == "E-GEOD-67297"),]
	#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-67297/?query=brown+adipose+tissue&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22

#get this dataset
E_GEOD_67297=ArrayExpress("E-GEOD-67297")

#see manufacturer
E_GEOD_67297@manufacturer

#As E_GEOD_67297 is an Affymetrix experiment, we can use RMA normalisation to process the raw data, please read the rma function help.
library("affy")
AEsetnorm = rma(E_GEOD_67297)

expr = assayDataElement(E_GEOD_67297,"G")

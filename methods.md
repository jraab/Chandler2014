**Computational Analysis:**
 
All analyses were performed using R custom scripts. Frequency of mutations in ARID1A and/or PIK3CA were identified by querying the [cBioPortal](http://www.cbioportal.org/public-portal/) on 2013-11-13 using the R API. Citation Cerami et al Cancer Discov. 2013 & Gao et al Sci. Signal. 2013. Total number of mutations for each tumor were calculated based on TCGA Somatic Mutation reports downloaded from the TCGA Data portal on 2013-10-29 for cancers with greater than 5% mutations in ARID1A and with public mutation data. These were stad-Stomach Adenocarcinoma, skcm - Skin Cutaneous Melanoma, luad - Lung Adenocarcinoma, lusc - Lung Squamous Cell Carcinoma, ucec - Uterine Corpus Endometrial Carcinoma. MSI data was downloaded from [Kim et al 2013](http://www.sciencedirect.com/science/article/pii/S0092867413012919) supplementary figure 1. TCGA barcodes in this figure for UCEC and the TCGA barcodes from the UCEC Somatic Mutation file referenced above were used to identify which samples contained mutations in ARID1A. All code used in these analysis are provided at http://github.com/jraab/Chandler2014.git


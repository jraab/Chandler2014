#quick look at mutation frequency in UCEC 
#data downloaded 2013-10-23 from TCGA hub as somatic mutations
library(plyr)
library(ggplot2)
setwd('/media/data/tcga/data/ucec/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/')
bi <- read.table('broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf', header=T, sep='\t')
dim(bi)
names(bi)
head(bi[,'Tumor_Sample_Barcode'])
arid1a_muts <- bi[bi$Hugo_Symbol %in% 'ARID1A', ]
arid1a_muts_id <- arid1a_muts$Tumor_Sample_Barcode
uniq_arid1a_id <- unique(arid1a_muts_id)
pik3ca_muts <- bi[bi$Hugo_Symbol %in% 'PIK3CA', ]
pik3ca_muts_id <- pik3ca_muts$Tumor_Sample_Barcode
uniq_pik3ca_id <- unique(pik3ca_muts_id)
uniq_total_id <- unique(bi$Tumor_Sample_Barcode)

totals_by_sample <- ddply(bi, .(Tumor_Sample_Barcode), nrow )
colnames(totals_by_sample) <- c('Tumor_Sample_Barcode', 'total_mutations' )
summary(totals_by_sample$total_mutations)
arid1a_samples <- totals_by_sample[totals_by_sample$Tumor_Sample_Barcode %in% uniq_arid1a_id, ]
pik3ca_samples <- totals_by_sample[totals_by_sample$Tumor_Sample_Barcode %in% uniq_pik3ca_id, ]
notarid1a_samples <- totals_by_sample[!(totals_by_sample$Tumor_Sample_Barcode %in% uniq_arid1a_id), ]

arid1a_only <- totals_by_sample[totals_by_sample$Tumor_Sample_Barcode %in% uniq_arid1a_id &
                                 !( totals_by_sample$Tumor_Sample_Barcode %in% uniq_pik3ca_id), ]
pik3ca_only <- totals_by_sample[totals_by_sample$Tumor_Sample_Barcode %in% uniq_pik3ca_id & 
                                  !(totals_by_sample$Tumor_Sample_Barcode %in% uniq_arid1a_id),]
both <- totals_by_sample[totals_by_sample$Tumor_Sample_Barcode %in% uniq_arid1a_id &
                           (totals_by_sample$Tumor_Sample_Barcode %in% uniq_pik3ca_id), ]
none <- totals_by_sample[!(totals_by_sample$Tumor_Sample_Barcode %in% uniq_arid1a_id) &
                           !(totals_by_sample$Tumor_Sample_Barcode %in% uniq_pik3ca_id), ]
df2 <- data.frame(mut_num=c(arid1a_only$total_mutations, 
                                pik3ca_only$total_mutations, 
                                both$total_mutations, 
                                none$total_mutations), 
                  group =c(rep('arid1a', nrow(arid1a_only)),
                               rep('pik3ca', nrow(pik3ca_only)), 
                               rep('both', nrow(both)), 
                               rep('none', nrow(none)) ))

ggplot(df2, aes(x=group, y=mut_num)) + geom_boxplot(outlier.size=0) + geom_jitter()
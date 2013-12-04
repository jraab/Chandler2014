#script to analyze data from 
#The Landscape of Microsatellite Instability in Colorectal and Endometrial Cancer Genomes
#Kim et al 2013
#see supplemental table 1 for source data. 
date()
library(plyr)
library(ggplot2)
library(gridExtra)
setwd('/media/data/tcga/data')
msi = read.csv('msi_published.csv', header=T) #this is supplemental table 1
bi <- read.table('ucec/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf', header=T, sep='\t')

getid <- function(x) { 
#split the tcga id and return the important part
#input is the three part ids from MSI sample
  parts = unlist(strsplit(as.character(x), "-"))
  id = paste(parts[2], parts[3], sep= "-")
  return(id) 
}

#parse out which samples are mutant and which are WT for ARID1A and fix the barcode ID
arid1a_muts <- bi[bi$Hugo_Symbol %in% 'ARID1A', ]
arid1a_muts_id <- arid1a_muts$Tumor_Sample_Barcode
uniq_arid1a_id <- unique(arid1a_muts_id)
uniq_total_id <- unique(bi[!bi$Tumor_Sample_Barcode %in% arid1a_muts_id, ]$Tumor_Sample_Barcode )
arid1a_df = data.frame(id = unlist(llply(uniq_arid1a_id, getid)), 
                       status = rep('mut', length(uniq_arid1a_id)))
wt_df = data.frame(id = unlist(llply(uniq_total_id, getid)), 
                   status=rep('wt', length(uniq_total_id)))
wt_df = wt_df[!(wt_df$id%in% uniq_arid1a_id),]
df = rbind(arid1a_df, wt_df)

# select the parts of the MSI frame I care about and fix the IDs
#msi_ucec = msi[msi$Histologic.type == 'Endometrioid',] # all endometroid tumors are ucec
msi_ucec = msi
msi_ucec = msi_ucec[!is.na(msi_ucec$Histologic.type), ]
msi_ucec_ids = as.character(msi_ucec$TCGA.identifier)
f = unlist(llply(msi_ucec_ids, getid))
msi_ucec$id = f
msi_ucec_info = msi_ucec[c('MSI.category', 'MLH1.silencing', 'MSI.events', 'Stage', 'Histologic.type', 'id')]

#where possible - assign wt/mut status to the samples int he MSI data set
final = merge(x=msi_ucec_info, y=df, by='id', all=FALSE)

#create some summary statistics by category
counts_by_status = ddply(final, .(status, MSI.category), summarize, count=length(MSI.category))
counts_by_status_norm = ddply(counts_by_status, .(MSI.category), transform, percent=count/sum(count))


#plotting
a <- ggplot(counts_by_status, aes(x=MSI.category, y=count, fill=status))
a <-a + geom_bar(stat='identity') + theme_classic() + theme(legend.position='none')
a <- a+ xlab('') + ylab('Number of Tumors') + labs(fill='')
b <- ggplot(counts_by_status_norm, aes(x=MSI.category, y=percent, fill=status)) + geom_bar(stat='identity')
b <- b + geom_bar(stat='identity') + theme_classic()
b <- b + xlab('') + ylab('Percent of Total MSI Samples') + labs(fill="")
pdf('/media/data/tcga/plots/ucec_msi_plots.pdf', height=5, width=10)
grid.arrange(a,b, ncol=2)
dev.off()

#note that ucec has two subtypes (plus the mixed category)
#again ARID1A not mutated in the Serous at an appreciable level )

pdf('../plots/ucec_subtypes_byMutation.pdf')
ggplot(final, aes(x=status, fill=Histologic.type)) + geom_bar(stat='bin') + theme_classic()
dev.off()
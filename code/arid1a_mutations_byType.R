#script to get mutation data for multiple cancer types. 
#check how many mutations per sample /plotted based on ARID1A status
library(plyr)
library(ggplot2)
setwd('/media/data/tcga/data/') # this needs adjusted to your data set
#setwd(~/PATH/TO/DATA/)
dir  <- '/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/'
#cancers <- c('blca', 'cesc', 'lgg', 'luad', 'lusc', 'skcm', 'stad', 'ucec')
#I'm using only a subset of these, becasue of issues with publicness of the data
cancers <- c('luad', 'lusc', 'skcm', 'stad', 'ucec')
df <- data.frame(cancer=NULL, arid1a_status=NULL, mutation_number=NULL)
pval <- data.frame(pval=rep(NA,length(cancers)))
rownames(pval) <- cancers
for (i in cancers) { 
  bi <- read.table(paste(i,dir,'/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf', sep='')
                         , header=T, sep='\t')
  arid1a_muts <- bi[bi$Hugo_Symbol %in% 'ARID1A', ]
  uniq_arid1a_id <- unique( arid1a_muts$Tumor_Sample_Barcode )
  uniq_total_id <- unique(bi$Tumor_Sample_Barcode)
  totals_by_sample <- ddply(bi, .(Tumor_Sample_Barcode), nrow )
  colnames(totals_by_sample) <- c('Tumor_Sample_Barcode', 'total_mutations' )
  arid1a_only <- totals_by_sample[totals_by_sample$Tumor_Sample_Barcode %in% uniq_arid1a_id, ] 
  none <- totals_by_sample[!(totals_by_sample$Tumor_Sample_Barcode %in% uniq_arid1a_id), ]
  df2 <- data.frame(cancer=rep(i, nrow(arid1a_only)+nrow(none)), 
                    arid1a_status=c(rep('mut', nrow(arid1a_only)), rep('wt',nrow(none) )), 
                    mutation_number=c( arid1a_only$total_mutations, none$total_mutations) )
  df <- rbind(df, df2 )    
  pval[i,'pval'] <- wilcox.test(arid1a_only$total_mutations, none$total_mutations)$p.value
}

print(xtable(x=pval ), digits=5, type='html', file='~/Dropbox/test.html')
print(ddply(df, .(cancer, arid1a_status), nrow))
medianDf <- ddply(df, .(cancer, arid1a_status), summarise, median=median(mutation_number))
medianDf <- medianDf[with(medianDf, order(medianDf$median, decreasing=T)), ]
s <-subset(medianDf, subset=medianDf$arid1a_status=='mut')
medianDf$cancer <- factor(as.character(medianDf$cancer), levels=unique(s$cancer))
medianDf$arid1a_status <- factor(as.character(medianDf$arid1a_status), levels=c('mut', 'wt'))
f <- ggplot(medianDf, aes(x=cancer, y=median, colour=arid1a_status))
f + geom_point(size=15, shape=95) + theme_classic() + xlab('') 

#this plot is ok, but can't do the jittering with the points easily. 
df$cancer <- factor(as.character(df$cancer), levels=unique(s$cancer)) #calc'd above
df$arid1a_status <- factor(as.character(df$arid1a_status), levels=c('mut', 'wt'))
g <- ggplot(df, aes(x=cancer, y=mutation_number, fill=arid1a_status)) + 
  geom_boxplot(outlier.size=NA, position=position_dodge(width=1))
g + scale_y_continuous(limits=c(0,4000)) + theme_classic() + xlab('') 



# this plot is the winner. 
pdf('/media/data/tcga/plots/cross_cancer_mutationrate_byArid1a_status.pdf')
h <- ggplot(df, aes(x=arid1a_status, y=mutation_number, fill=arid1a_status)) + 
geom_boxplot(alpha=0.7, outlier.size=NA) + facet_grid(~cancer) 
h <- h + geom_jitter(alpha=0.4, aes(colour=arid1a_status)) + scale_y_continuous(limits=c(0,4000))
h + theme_bw() + xlab('') + ylab('Total Mutations Per Sample') + 
  theme(legend.position='none', strip.background=element_rect(fill='white'), panel.grid=element_blank())
dev.off()
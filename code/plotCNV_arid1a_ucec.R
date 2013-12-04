#what to do with copy number
setwd('/media/data/tcga/data/ucec/CNV_SNP_Array/BI__Genome_Wide_SNP_6//Level_3')
meta <- read.table('/media/data/tcga/data/ucec/METADATA/BI__Genome_Wide_SNP_6/broad.mit.edu_UCEC.Genome_Wide_SNP_6.sdrf.txt', 
                   header=T, sep='\t')
mapping <- data.frame(barcode = meta$Comment..TCGA.Barcode., 
                      fn = paste(meta$Hybridization.Name, '.hg19.seg.txt', sep=''))
df <- data.frame(barcode=NA, chr=NA, start=0, end=0,  seg_mean=0, sample=NA)

dir  <- '/media/data/tcga/data/ucec/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/'
bi <- read.table(paste(dir,'/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf',sep='')
                 , header=T, sep='\t')
arid1a_muts <- bi[bi$Hugo_Symbol %in% 'ARID1A', ]
arid1a_mutants <- unique( arid1a_muts$Tumor_Sample_Barcode )

getCode <- function(x) {
  z <- strsplit(as.character(x), '-')[[1]]
  id <- paste(z[2:3], collapse='-')
  return(id)
}
arid1a_mutants <- unlist(llply(arid1a_mutants, getCode))


for (i in 1:nrow(mapping)) { 
  row <- mapping[i,]
  name <- as.character(row[1,2])
  possibleError <- tryCatch(read.table(name), error=function(e) e)
  if (inherits(possibleError, 'error')) next
  fn <- read.table(name, header=T)
  tmp<- strsplit(as.character(row[1,1]), split='-')[[1]]
  id <- getCode(row[1,1])
  code <- tmp[4]
  if (code >= 10 ) {
    sample <- 'normal'
  }
  else {
    if (id %in% arid1a_mutants) {
      sample <- 'tumor_mut'
    }
    else {
      sample <- 'tumor_wt'
    }
  }
  df2 <- data.frame(barcode=as.character(row[1,1]),
                    chr=fn$Chromosome, 
                    start=fn$Start, 
                    end=fn$End, 
                    seg_mean=fn$Segment_Mean, 
                    sample=sample )
  df <- rbind(df, df2) 
}

df <- df[!is.na(df$chr),]

g<- ggplot(df, aes(x=start,y=seg_mean,colour=sample)) 
g<- g + facet_grid(sample~chr,scales='free') + theme(panel.margin=unit(0,'lines'))
g+ geom_point(size=0.75, alpha=0.25) + coord_fixed(ratio=1/10)
# this looks ok, but I don't see evidence that ARID1A mutant tumors are any different
#part of the problem is that I only get a point if there is a change, so when I have tons of samples
# in one group and few in the other it is difficult to compare. 

#what if I look at just how many copy number changes are in each group
d <- ddply(df, .(barcode), summarise, count=length(seg_mean), sample=unique(sample))
ggplot(d, aes(x=sample, y= count, fill=sample)) + geom_boxplot(outlier.size=NA, notch=TRUE) + scale_y_continuous(limits=c(0,2000))
wilcox.test(subset(d, subset=d$sample =='normal')$count, subset(d, subset=d$sample=='tumor_wt')$count)
wilcox.test(subset(d, subset=d$sample =='normal')$count, subset(d, subset=d$sample=='tumor_mut')$count)
wilcox.test(subset(d, subset=d$sample =='tumor_mut')$count, subset(d, subset=d$sample=='tumor_wt')$count)
#although both tumors samples have more changed copy number segments, there is no difference between arid1a mut or w
#what if I look at only copy number changes greater or less than some amount
#abs(seg_mean >0.5)
d <- ddply(df, .(barcode), summarise, count=length(seg_mean[abs(seg_mean) > 0.5]), sample=unique(sample))
ggplot(d, aes(x=sample, y= count, fill=sample)) + geom_boxplot(outlier.size=NA, notch=TRUE) + scale_y_continuous(limits=c(0,300))
wilcox.test(subset(d, subset=d$sample =='normal')$count, subset(d, subset=d$sample=='tumor_wt')$count)
wilcox.test(subset(d, subset=d$sample =='normal')$count, subset(d, subset=d$sample=='tumor_mut')$count)
wilcox.test(subset(d, subset=d$sample =='tumor_mut')$count, subset(d, subset=d$sample=='tumor_wt')$count)

#abs(seg_mean > 1)
d <- ddply(df, .(barcode), summarise, count=length(seg_mean[abs(seg_mean) > 1]), sample=unique(sample))
ggplot(d, aes(x=sample, y= count, fill=sample)) + geom_boxplot(outlier.size=NA, notch=TRUE) + scale_y_continuous(limits=c(0,300))
wilcox.test(subset(d, subset=d$sample =='normal')$count, subset(d, subset=d$sample=='tumor_wt')$count)
wilcox.test(subset(d, subset=d$sample =='normal')$count, subset(d, subset=d$sample=='tumor_mut')$count)
wilcox.test(subset(d, subset=d$sample =='tumor_mut')$count, subset(d, subset=d$sample=='tumor_wt')$count)

#performs slightly worse at abs(2)
d <- ddply(df, .(barcode), summarise, count=length(seg_mean[abs(seg_mean) > 2]), sample=unique(sample))
ggplot(d, aes(x=sample, y= count, fill=sample)) + geom_boxplot(outlier.size=NA, notch=TRUE) + scale_y_continuous(limits=c(0,300))
wilcox.test(subset(d, subset=d$sample =='normal')$count, subset(d, subset=d$sample=='tumor_wt')$count)
wilcox.test(subset(d, subset=d$sample =='normal')$count, subset(d, subset=d$sample=='tumor_mut')$count)
wilcox.test(subset(d, subset=d$sample =='tumor_mut')$count, subset(d, subset=d$sample=='tumor_wt')$count)


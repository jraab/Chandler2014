### This script is to generate plots for Ron's PIK3CA/ARID1A tumor experiments
### Most of the data comes from cBio portal, but a separate figure should be generated to deal 
### with data from the ovarian clear cell carcinoma science paper. I collated this data from 
### the supplementary data from that paper. 
#initial work for this scrip is in the same directory analysis_initial.Rmd

date()
library(plyr)
library(reshape2)
library(ggplot2)
library(cgdsr)
library(gridExtra)

setwd('~/Dropbox/github/labwork/mutations/')


###function defs
replaceNaN <- function(cell){ 
  if (is.na(cell) ) {return(NA)}
  if (cell == 'NaN' ) {cell=NA}
  return(cell)
}

mutType <- function(x) { 
  if (x[1] == 1 & x[2] ==0){ a <- 'ARID1A'}
  if (x[1] == 0 & x[2] ==1) {a <- 'PIK3CA'}
  if (x[1] ==1 & x[2] ==1) { a <- 'BOTH'}
  if (x[1] ==0 & x[2] ==0 ) {a <- 'NONE'}
  return(a)
}

find_dups <- function(list_of_ids) { 
  index <- list()
  for (i in list_of_ids) { 
    p1 <- strsplit(i, '_')[[1]][1] 
    if (length(grep(p1, list_of_ids)) >1 ){ 
      patt <- 'pub'
      id <- as.character(grep(patt,list_of_ids,value=TRUE, invert=TRUE))
      index <- c(index,id)
    }
    else {
      i <- as.character(i)
      index <- c(index,i)
    }
  }
  return(unique(unlist(index)))
}
#CC3300'
values=c( '#99CCFF','#666666','#CCCCCC')
########################
#occc graph
pik3ca_scisupp <- read.table(file='~/Dropbox/github/labwork/mutations/data/PIK3CA_ARID1A_OCCC_Science_supp.csv', header=T, sep='\t', stringsAsFactors=FALSE)
pik3ca_scisupp$ARID1A.Mutation <- gsub(" ", "", pik3ca_scisupp$ARID1A.Mutation)
pik3ca_scisupp$PIK3CA.Mutation <- gsub(" ", "", pik3ca_scisupp$PIK3CA.Mutation)
both <- ifelse(pik3ca_scisupp$ARID1A.Mutation=='Y' & pik3ca_scisupp$PIK3CA.Mutation=='Y', 1, 0)
arid1a_occc <- ifelse(pik3ca_scisupp$ARID1A.Mutation=="Y", 1, 0 )
pik3ca_occc <- ifelse(pik3ca_scisupp$PIK3CA.Mutation=="Y", 1, 0 )
df <- data.frame(arid1a=arid1a_occc, pik3ca=pik3ca_occc)
df2 <- apply(df, MARGIN=1, FUN=mutType)
df3 <- as.data.frame(table(df2))
total <- sum(df3$Freq)
df3 <- df3[!df3$df2 == 'NONE',]
df3$Freq <- df3$Freq/total *100
df3$df2 <- factor(df3$df2, levels=as.character(c('BOTH', 'ARID1A', 'PIK3CA')))
df3 <- df3[with(df3, order(df3$df2)), ]
occc_data <- df3
#cancer <- 'Ovarian Clear Cell Carcinoma - Jones 2010'
cancer <- 'OCCC'
g <- ggplot(df3, aes(x=cancer, fill=df2, y = Freq)) + geom_bar(stat='identity')
g <- g + xlab('') + ylab('Percent of Tumors with Mutations')+ggtitle('Mutation frequencey in OCCC')
g <- g + coord_flip(ylim=c(0,100)) + labs(fill='')
g <- g + theme(axis.ticks.x=element_blank(), axis.text.x = element_blank())

occc_graph <- g + theme_bw() + scale_fill_manual(values=values) + coord_fixed(ratio=.1/1, ylim=c(0,80))  + 
  theme(legend.position='none') 
#################
#data for IHC ARID1A and PIK3CA sequencing - numbers from RC. 

#total <- 42
#arid1a_def_pik3_mut <- 12/total * 100
#arid1a_def <- 11/total * 100 # not including the 12 above
#pik3_mut <- 5/total * 100
#normal <- 14/total* 100
#expt <- 'OCCC IHC Data - NEJM'
#df <- data.frame(cancer=expt, Both = arid1a_def_pik3_mut, ARID1A = arid1a_def, PIK3CA = pik3_mut)
#ihcdata <- df
#df_m <- melt(df, id.vars='cancer')
#g <- ggplot(df_m, aes(x=expt, y=value, fill=variable)) + geom_bar(stat='identity') + coord_fixed(0.1/1, ylim=c(0,80))
#g <- g + theme_bw() + scale_fill_manual(values=values) 
#ihc_graph <- g + ylab('') + xlab('') + ggtitle('IHC data') + theme(legend.position='none')
#####################
#combine occc and ihc data ito one panel. 
#ihcrow <- c(as.character(ihcdata$cancer), ihcdata$Both, ihcdata$ARID1A, ihcdata$PIK3CA)
occcrow <- c(cancer, occc_data$Freq[1] , occc_data$Freq[2], occc_data$Freq[3])
df <- data.frame(matrix(occcrow,nrow=1) )
#print(df)
colnames(df) <- c('Cancer', 'BOTH', 'ARID1A', 'PIK3CA')
df_m <- melt(df, id.vars='Cancer')
df_m$value <- as.numeric(df_m$value)
#g <- ggplot(df_m, aes(x=Cancer, y = value, fill=variable)) + geom_bar(stat='identity')
#g <- g +scale_fill_manual(values=values) + theme_bw() + theme(axis.text.x=element_text(angle=90)) 
#g <- g + xlab('') + ylab('') + ylim(c(0,80)) + theme(legend.position='none') 
#combined <- g+ ggtitle('Ovarian Clear Cell Mutations') + coord_flip() + ylim(c(0,80)) + theme(panel.border=element_blank())
#combined <- combined + theme(axis.line=element_line(colour='black'))

############
library(ggplot2)
library(reshape2)
library(plyr)
library(cgdsr)

replaceNaN <- function(cell){ 
  if (is.na(cell) ) {return(NA)}
  if (cell == 'NaN' ) {cell=NA}
  return(cell)
}

find_dups <- function(list_of_ids) { 
  index <- list()
  for (i in list_of_ids) { 
    p1 <- strsplit(i, '_')[[1]][1] 
    if (length(grep(p1, list_of_ids)) >1 ){ 
      patt <- 'pub'
      id <- as.character(grep(patt,list_of_ids,value=TRUE, invert=TRUE))
      index <- c(index,id)
    }
    else {
      i <- as.character(i)
      index <- c(index,i)
    }
  }
  return(unique(unlist(index)))
}

getInfo <- function(study_id, genes) { 
  study_id <- study_id
  caseList = paste(study_id, '_cnaseq', sep='') #cases with gistic and mutation data
  study_info <- getCaseLists(mycgds, study_id)[1] 
  if (caseList %in% unlist(study_info)){
    gistic <- getProfileData(mycgds, genes=genes, geneticProfiles = paste(study_id, '_gistic', sep=''), caseList=caseList)
    mutation <- getProfileData(mycgds, genes=genes, geneticProfiles = paste(study_id, '_mutations', sep=''), caseList=caseList)
    m <- merge(gistic, mutation, by = 'row.names')
    return(study_id=m)
  }
  else{ 
    string <- paste('Failed-NoValid _cnaseq set', study_id, sep=' ' )
    print(string)
  }
}  

isMut <- function(row) {  
  row <- as.list(row) # no clue why this is necessary, but otherwise can't get the numeric line to take
  row[2] <- as.numeric(row[2])
  mut <- ifelse(row[2] < -1 | row[2] > 1 | !is.na(row[3]), 1, 0 )
  return(mut)
}

mutType <- function(x) { 
  if (x[1] == 1 & x[2] ==0){ a <- 'ARID1A'}
  if (x[1] == 0 & x[2] ==1) {a <- 'PIK3CA'}
  if (x[1] ==1 & x[2] ==1) { a <- 'BOTH'}
  if (x[1] ==0 & x[2] ==0 ) {a <- 'NONE'}
  return(a)
}

mapNames <- function(study_id) { 
  #maps study_id to printable names
  require(cgdsr)
  mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
  studyMapping <- getCancerStudies(mycgds)[,c(1:2)]
  name <- studyMapping[studyMapping$cancer_study_id %in% study_id, ]$name
  name <- unlist(strsplit(name, "\\(", ))[[1]]
  return(name)
}

mycgds <- mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
studies <- getCancerStudies(mycgds)[,c(1,2)] 
studies_noDup <- find_dups(studies$cancer_study_id)
studies <- studies[ studies$cancer_study_id %in% studies_noDup, ]
studies <- studies[grep(pattern='ccle', x=studies$cancer_study_id, invert=T),] #cell lines aren't important and screw things 
#b/c there is not GISTIC data. 
genes <- c('ARID1A', 'PIK3CA')
#foreach study - get the total number of mutations, #ARID1A, #PIK3CA, #BOTH
#getInfo returns a dataframe where each row is a tumor with mutation info
out <- data.frame(matrix(NA, nrow=length(studies[,1]), ncol=5))
pval <-data.frame(study=as.factor(studies[,1]), pval=rep(NA), or=rep(NA))
colnames(out) <- c('Cancer', 'BOTH', 'ARID1A', 'PIK3CA', 'NONE' )
j <- 1
for (i in studies[,1]){
  cases <- getCaseLists(mycgds, cancerStudy=i)
  cases <-getInfo(i,'ARID1A') #dummy df so that i can get number of rows
  if (is.numeric(nrow(cases))) {
    bothMuts <- data.frame(matrix(nrow=nrow(cases)))
  }
  else{
  next
  }
  for (a in genes){ 
    df <- getInfo(i, a)
    muts <- apply(df, 1, isMut)
    bothMuts[,eval(a)] <- muts
    
  }
  bothMuts$cancer <- rep(eval(i))
  bothMuts<- bothMuts[,-1]
  studyName <- bothMuts[1,3]
  bothMuts <- bothMuts[,-3]
  # bothMuts <- data.frame(apply(bothMuts, c(1,2), replaceNaN))
  df2 <- apply(bothMuts, 1, mutType)
  df3 <- as.data.frame(table(factor(df2, levels=c('BOTH', 'ARID1A', 'PIK3CA', 'NONE'))))
  arid1a_num <- df3[df3$Var =='ARID1A', ]
  pik3ca_num <- df3[df3$Var =='PIK3CA', ]
  both_num <- df3[df3$Var =='BOTH', ]
  none_num <- df3[df3$Var =='NONE', ]
  m <- matrix(c(both_num$Freq, arid1a_num$Freq, pik3ca_num$Freq, none_num$Freq), nrow=2)
  ftest <- fisher.test(m)
  p <- ftest$p.value
  or <- ftest$estimate
  pval[pval$study==i,]['pval'] <- p
  pval[pval$study==i, ]['or'] <- or
  total <- sum(df3$Freq)
  df3$Freq <- df3$Freq/total *100
  #out[j ,] <- c(studyName , as.numeric(df3$Freq))
  out[j, ] <- c(toupper(studyName), as.numeric(df3$Freq)) 
  j <- j +1
}

out <- out[complete.cases(out),]
print(out)

out <- out[,1:4] # remove the none category
out[,2] <- as.numeric(out[,2])
out[,3] <- as.numeric(out[,3])
out[,4] <- as.numeric(out[,4])
out2 <- out[with(out, order(out[,'BOTH'], decreasing=T) ), ]
#out2$Cancer <- llply(out2$Cancer, .fun=mapNames) #set this to get the pretty name for study
#pval$study <- llply(pval$study, .fun=mapNames) #set this to get the pretty name for study
pval <- pval[!is.na(pval$pval),]
pval$study <- sapply(pval$study, as.character)
write.table(x=pval, file='pvals.tsv', sep='\t', row.names=FALSE)
##########need the following to 
#plot as all cancers on 1 graph
#cancer <- 'Ovarian Clear Cell Carcinoma (Jones ,Science)'
#occc_row <- c(cancer, as.numeric(df3$Freq[1]), as.numeric(df3$Freq[2]), as.numeric(df3$Freq[3]))
#t <- rbind(out, occc_row)
#########################

out2$Cancer <- factor(out2$Cancer, levels=unique(as.character(out2$Cancer)))
out2$Cancer <- gsub("_TCGA", "", out$Cancer)#needed to clean up the label
out_m <- melt(out2, id.vars='Cancer')
g <- ggplot(out_m, aes(x=Cancer, y=value, fill=variable)) + geom_bar(stat='identity')
g <- g + theme_bw() + xlab('') +ylab('% of Tumors Mutated')  + labs(fill="")
g <- g  + theme(legend.position=c(1,0),  legend.justification=c(1,0), legend.background=element_blank()) 
cbio_graph <- g + ggtitle('Mutation Frequency From cBio') + coord_flip(ylim=c(0,80) ) + scale_fill_manual(values=values)
cbio_graph <- cbio_graph + theme(axis.line=element_line(colour='black'), panel.border=element_blank())
#ggsave('occ_graph.pdf', occc_graph, height=10, width=5, path='~/Dropbox/github/labwork/mutations/plots/')
#ggsave('cbio_graph.pdf', cbio_graph, width=10, path='~/Dropbox/github/labwork/mutations/plots/')


#out2$Cancer <- factor(out2$Cancer, levels=unique(as.character(out2$Cancer)))
#out_m <- melt(out2, id.vars='Cancer')
#g <- ggplot(out_m, aes(x=Cancer, y=value, fill=variable)) + geom_bar(stat='identity')
#g <- g + theme_bw() + coord_flip() + xlab('') +ylab('% of Tumors Mutated') + scale_y_continuous()


write.table(file='/media/data/tcga/Chandler2014/data/arid1a_pik3ca_overlaps.tsv', x=out2, sep='\t', col.names=TRUE, row.names=FALSE)

##################
#final plots using fake facets
#setwd('~/Dropbox/github/labwork/mutations/plots/')
#pdf('final.pdf', width=8.5, height=11,paper='letter' )
#final <- grid.arrange(cbio_graph, combined, nrow=2, heights = c(3, 1))
#final
#dev.off()

out_m$facet <- 'facet2'
df_m$facet <- 'facet1'
r <- rbind(df_m, out_m )
g <- ggplot(r, aes(x=Cancer, y=value, fill=variable)) + facet_grid(.~facet, scale='free_x', space='free_x')
g <- g + geom_bar(position='stack', stat='identity') + scale_fill_manual(values=values) + theme_bw() 
g <- g + theme(axis.text.x=element_text(angle=90, hjust=1), 
               axis.text.y=element_text(angle=90),
               strip.text.x=element_blank(), 
               strip.background=element_blank(), 
               legend.position=c(.9,.9),
               legend.background=element_blank(),
               axis.line=element_line(colour='black'), 
               panel.border=element_blank(), 
               panel.margin = unit(2, "lines"),
               panel.grid = element_blank()
               ) 

g <- g + ylab('Percent Altered') +labs(fill='') + xlab('') 
ggsave('/media/data/tcga/Chandler2014/plots/arid1a_pik3ca_overlap_crossCancer.pdf', g, height=6, width=6, units='in')

 

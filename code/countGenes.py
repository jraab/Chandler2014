#!/usr/bin/python

# script to scan through all files in an directory and get information about expression of genes of interest.
# if I've done this right it should be able to go through a tcga directory and output a pandas dataframe of sample X gene expression values
import os
import re
import pandas as pd
import numpy as np
patt = re.compile('gene_normalized')
seqdir ='/magnuson-lab/jraab/tcga/data/hcc/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/' 
files = [f for f in os.listdir(seqdir) if 'gene_normalized' in f ]  
expression = {}
genes= ['ARID1A', 'ARID1B', "ARID2", 'SMARCA4', 'SMARCA2', 'CTNNB1', 'PIK3CA'] 
barcode =set() 
for name in files:  
   barcode.add(name.split('__')[2])

df = pd.DataFrame( index = barcode, columns = genes)
for f in files:
   with open(seqdir+f) as fn:  
      for l in fn:  
         for g in genes: 
            g_patt = re.compile(r"^"+g)  
            if g_patt.match(l.split('\t')[1]):  
               barcode, gene, expr = l.split('\t')
               gene = gene.split('|')[0] 
               df[gene][barcode] = expr 


df.to_csv('test.txt')

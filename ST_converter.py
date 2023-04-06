# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 08:56:24 2019

@author: mabro

script to reformat dataset to the suppl table format

"""
import os
import sys
import pandas as pd
import numpy as np
sys.path.append(os.path.abspath('Z:/mabro/QTs/paper/scripts/QT-GWAS/utilities')) #path to  where your config file is located containing all path variables (PT_<yourfile>)
from config_suppltable_conv import *

data = pd.read_csv(PT_data,sep = '\t', header = None)
data.columns = ['index','index2','association', 'chr', 'traitnr', 'trait', 'trait2', 
                     'trait3', 'start', 'stop', 'folder', 'snp', 'pval', 'odds' ,
                     'kinship_AV', 'kinship_SD', 'conttable','filled_pval', 'chrom', 'loc', 
                     'filled_loc', 'filled_chrom', 'filled_trait']

loci = data[data['folder'].notnull()].index.tolist() #get indices where locus information is stored
tab = data[['trait3', 'snp', 'pval','start','stop']] # these columns contain the data we need for the suppl table
tab['genes'] = ''


df_subset = data[1:38]

#function to retrieve genes from  a certain association in 1 string instead of splitted over multiple lines
def getgenes(df_subset):
    genes = df_subset[df_subset['chr'].notnull()] #in the filtered fulldataset the genes involved in a certain association 
    #are currently stored under the column 'chr', select all non empty values to get all the rows containing genes involved in the current association
    if len(genes) == 0:
        return('')
    else:
        genes_spl = genes['chr'].str.split(' ', expand = True)    
        gene_list = genes_spl[0].to_list()
        s = ''
        for gene in gene_list: 
            s += gene + ', ' #make a string of all the genes in the current association
        s = s.rstrip(', ') #remove last comma
        return (s)

#for all loci retrieve the associated genes as a string and fill this string into the genes column
for i in range (len(loci)):
    start = loci[i]+1 #index at which the subset starts: index first row below the row containing the current locus information)
    if i != len(loci)-1: #if not at last locuss
        stop = loci[i+1]-1 #index at which the subset stops: index of last row before a new locus row 
    else: #if at last locus
        stop = len(data)-1 #stop is the end of the data frame
    subset = data[start:stop]
    s = getgenes(subset)
    tab.iat[loci[i],5] = s #fill in the string containing all associated genes into column 5

#select only the rows containing locus information    
tab = tab.loc[loci]
tab.columns = ['trait', 'SNP', 'pval','start','stop', 'associated genes']
tab.to_csv(PT_output,sep = '\t', index = False) #store into a file

# -*- coding: utf-8 -*-
'''
Created on Thu Aug 22 11:27:27 2019

@author: mabro

fulldataset analyses
'''
#import modules
import os
import sys
import pandas as pd
import numpy as np
sys.path.append(os.path.abspath('Z:/mabro/QTs/paper/scripts/QT-GWAS/utilities')) #path to  where your config file is located containing all path variables (PT_<yourfile>)
from config_fulldata_analys import *


#read in data
fulldata = pd.read_csv(PT_data,sep='\t',header = None)
fulldata.columns = [['index','association', 'chr', 'traitnr', 'trait', 'trait2', 
                     'trait3', 'start', 'stop', 'folder', 'snp', 'pval', 'odds' ,
                     'kinship_AV', 'kinship_SD', 'conttable']] #assign column names

loci_indices = fulldata[fulldata['folder'].notnull()].index.tolist() #check at which lines we have locus information based on the 11th column: if there is value here, we are at the position of a new locus followed by several lines with info on the associated genes
len(loci_indices) #131193 associations


# function to fill in values of a certain column in the dataframe
def fillin(df,col): 
    new = [df[col][0]] #start with the first locus pval
    for i in range(1,len(df.index)): #for all rows after the first locus row
        if pd.isnull(df[col][i]): #if there is no value take the value of the row above
            new.append(new[i-1])
        elif df[col][i] == '': #for strings empty cells are represented by empty strings ''
            new.append(new[i-1])
        else:
            new.append(df[col][i]) #if there is a value we are at a new locus row: take this value for the next x rows belonging to the assoc
    return(new) #return as a list so we can assign it to a new column
fulldata['filled_pval'] = fillin(fulldata,'pval') #we need to fill in the pval of the association for all rows belonging to the same locus


#make a dictionary with for each trait the chromosome, location and pvalue of the associated SNPs
d = {}
for index in loci_indices:
    trait = fulldata.iloc[index][6]
    snp = fulldata.iloc[index][10]
    spl = snp.split("_")
    chrom = int(spl[0])
    pos = int(spl[1])
    pval = fulldata.iloc[index]['filled_pval']
    if trait in d.keys():
        d[trait].append([chrom,pos,pval])
    else:
        d.update({trait:[[chrom,pos,pval]]})
len(d) #1957 features

def redundant(d): #remove all redundant snps 
    for key in d:
        snplist = d[key]
        backup = tuple(snplist)
        for i,element in enumerate(backup): #compare each couple of snp and pval in the list to the other snps
            chrom = element[0]
            loc = element[1]
            pval1 = element[2]
            for j,element2 in enumerate(backup):
                if j != i: #do not compare snp to itself
                    chrom2 = element2[0]
                    loc2 = element2[1]
                    pval2 = element2[2]
                    if chrom2 == chrom and abs(int(loc2)-int(loc)) < 20000: #chrom is same and snps are within 20kb --> select the lowest pval
                        #take snp with lowest pval
                        if float(pval1) < float(pval2): #retain the smallest pvalue
                            if element2 in snplist:
                                snplist.remove(element2)

                        elif float(pval2) < float(pval1):
                            if element in snplist:
                                snplist.remove(element)
    return(d)
test = redundant(d) #run the redundant function on the dictionary containing all SNPs and p-values per trait, it returns a dictionary containing only nonredundant SNPs
data = fulldata.copy()
spl = data[10].str.split('_',expand = True)   #split the column containing the SNP information into chromosome and location of the SNPs 
data['chrom'] = spl[0].astype(np.float64)  #convert strings to floats
data['loc'] = spl[1].astype(np.float64) #convert strings to floats
 

#fill in the chromosome, location and trait values
data['filled_chrom'] = fillin(data,'chrom')
data['filled_loc'] = fillin(data,'loc')
data['filled_trait'] = fillin(data,'trait3')

#select all rows that correspond to a nonredundant association and store into a list of smaller_dfs
smaller_dfs = []
for key in test.keys():
    snps = test[key]
    for snp in snps:
        sel = data.loc[(data['filled_trait'] == key) & (data['filled_chrom'] == snp[0]) & (data['filled_loc'] == snp[1]) & (data['pval'] == snp[2])]
        smaller_dfs.append(sel)

large_df = pd.concat(smaller_dfs, ignore_index=True) #concatenate the small dataframes into a big one: this is the filtered dataset
large_df.to_csv(PT_output,sep = '\t', header = None) #write results to csv
loci_indices = large_df[large_df['folder'].notnull()].index.tolist() #check at which lines we have locus information based on the 11th column: if there is value here, we are at the position of a new locus followed by several lines with info on the associated genes     
len(loci_indices) #45382 loci




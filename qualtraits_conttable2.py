'''
script to construct contigencytables for qualitative trait analyses:
for each trait and associated snp determine how many ecotypes have the snp and trait, how many ecotypes lack snp and trait and how many
have snp but no trait and vice versa
note: the rows in the snpfile (contains the 250K SNPs data - 1 = present in ecotype, 0 = absent) and the columns trtfile 
(contains the 6790 LC-MS features, 1 = present in ecotype, 0 = absent) represent the 183 ecotypes and have to be in the same order 
'''

#import modules
import pandas as pd
import numpy as np
import sys 
import os
sys.path.append(os.path.abspath('Z:/mabro/QTs/paper/scripts/QT-GWAS/utilities')) #path to  where your config file is located containing all path variables (PT_<yourfile>)
from config_conttable2 import *

snpfile = pd.read_csv(PT_snpfile, sep = '\t')
trtfile = pd.read_csv(PT_trtfile, sep = '\t')
fulldataset =  pd.read_csv(PT_fulldata, sep = '\t', header = None)
fulldataset.columns = ['Trait', 'SNP', 'p-value', 'odds', 'Kinship_Average', 'Kinship_SD']

def conttable(snpfile,snp,trtfile,trait):
    snpdata = snpfile[snp] #find column at which the snp is
    arr = np.array(trtfile) #convert to numpy array
    row = list(arr[:,0]).index(trait) # slice out first col --> contains the traits and find the row at which current trait is
    searchlist_trait = list(arr[row,1:])# count in this row how many ones
    zero_zero = 0 # zero in snp and zero in trait val
    zero_one = 0 #zero in snp and 1 in trait
    one_zero = 0 #1 in snp and 0 in trait
    one_one = 0 #1 in both
    for i,snpval in enumerate(snpdata):
        traitval = searchlist_trait[i]
        if snpval == 0 and traitval == 0: #ecotype does not have the snp nor the trait
            zero_zero += 1
        elif snpval == 0 and traitval == 1: #ecotype does not have the snp but does have the trait
            zero_one += 1
        elif snpval == 1 and traitval == 0: #ecotype has the snp but lacks the trait
            one_zero += 1
        elif snpval == 1 and traitval == 1: #ecotype has both snp and trait. 
            one_one += 1
    conttab = [[zero_zero,one_zero],[zero_one,one_one]]
    return(conttab)


fulldataset['conttable'] = fulldataset.apply(lambda x: conttable(snpfile = snpfile, snp = x['SNP'], trtfile = trtfile, trait = x['Trait']), axis=1)

fulldataset.to_csv(PT_output, sep = '\t')





    
    
    
    
    
    
    
    
    




        
    



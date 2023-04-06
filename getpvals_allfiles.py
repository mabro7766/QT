
'''
script to merge the the output from the analyze_regions4.pl containing gene info with the "conttable2_fulldataset.txt" that contains snps, 
pvalues, odds, kinships and conttingency tables per association: these values will be filled in in the "fulldata_genes.txt" dataset
'''

#load required modules
import pandas as pd
import sys 
import os
sys.path.append(os.path.abspath('Z:/mabro/QTs/paper/scripts/QT-GWAS/utilities')) #path to  where your config file is located containing all path variables (PT_<yourfile>)
from config_getpvals import *

# read in both datasets
data = pd.read_csv(PT_data, sep = '\t')
genes = pd.read_csv(PT_genes, sep = '\t') 
genes.columns = genes.columns.str.strip() #remove whitespaces from column names for easy access

#split the snp information into chromosome and location values
data[['chromosome', 'location']] = data['SNP'].str.split('_', expand = True) #split the SNP column into chromosome and location
data['location'] = data['location'].astype('int64') #convert location values to integers


# create all columns to fill in in the genes dataset
genes['SNP'] = ''
genes['p-value'] = ''
genes['odds'] = ''
genes['Kinship_Average'] = ''
genes['Kinship_SD'] = ''
genes['conttable'] = ''


inds = genes.loc[genes['folder'] == 'folder'].index.tolist()
for i in inds:
    trait = genes.iloc[i]['Trait']
    chromosome = genes.iloc[i]['chromosome']
    start = genes.iloc[i]['start']
    match = data.loc[(data['Trait'] == trait) & 
                     (data['chromosome'] == chromosome) & (data['location'] == start + 10000)].index.to_list()[0]
    
    genes.iat[i,9] = data.iloc[match]['SNP'] #fill in the snp value in the genes dataset in the correct row (given by i) and column (9th column)
    genes.iat[i, 10] = data.iloc[match]['p-value'] #fill in the p-value in the genes dataset in the correct row (given by i) and column (10th column)
    genes.iat[i, 11] = data.iloc[match]['odds'] #fill in the odds value in the genes dataset in the correct row (given by i) and column (11th column)
    genes.iat[i, 12] = data.iloc[match]['Kinship_Average'] #fill in the kinship_average value in the genes dataset in the correct row (given by i) and column (12th column)
    genes.iat[i, 13] = data.iloc[match]['Kinship_SD'] #fill in the kinship_sd value in the genes dataset in the correct row (given by i) and column (13th column)
    genes.iat[i, 14] = data.iloc[match]['conttable']

genes.to_csv(PT_output, sep = '\t')


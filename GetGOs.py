# -*- coding: utf-8 -*-
"""
script to add GOterm info to suppl table 2 QT-GWAS
Created on Thu May  6 13:05:27 2021

@author: mabro
"""

import os
import sys
import pandas as pd
import numpy as np
sys.path.append(os.path.abspath('Z:/mabro/QTs/paper/scripts/QT-GWAS/utilities')) #path to  where your config file is located containing all path variables (PT_<yourfile>)
from config_GetGOs import *
from itertools import chain


data = pd.read_csv(PT_data, sep = '\t')
data['genes'] = data['genes'].fillna('')


#first split the column containing genes into multiple rows, retaining all other info

def chainer(s):
    return list(chain.from_iterable(s.str.split(',')))

# calculate lengths of splits
lens = data['genes'].str.split(',').map(len)

# create new dataframe, repeating or chaining as appropriate
cols = data.columns.tolist()

d = {}
for col in cols:
    if col != 'genes':
        d.update({col:np.repeat(data[col], lens)})
    else:
        d.update({col: chainer(data[col])})                 
res = pd.DataFrame(d)
#get all the genes in the dataset to get GOterms from tair
genes = res.genes
#genes.to_csv(PT_genes, header = True) 


#make a column with all go terms for each gene, then again split this like above

GO = pd.read_csv(PT_GO, sep = ';') #data obtained from TAIR: contains GOterms for each gene in the QT-GWAS dataset, separated by ';’ 

#add GOslims category for each GOterm


def getGOs(gene):    
    store = []
    goterms = ''
    ind = GO.index[GO['Locus'] == gene].tolist()
    for i in ind:
        ID = str(GO.iloc[i]['GO ID'])
        store.append(ID)
    store = list(set(store)) #remove duplicate goterms 
    for element in store:
        goterms += element + ', '
    goterms = goterms[:-2] #remove last ','
    return (goterms)


res['genes'] = res['genes'].str.strip() #remove all whitespace
res['GOterms'] = res['genes'].map(lambda x:getGOs(x)) #get all GOterms for each gene and store in  'GOterms' column

#split the GOterms as we did for genes
lens = res['GOterms'].str.split(',').map(len)

cols = res.columns.tolist()
d = {}
for col in cols:
    if col != 'GOterms':
        d.update({col:np.repeat(res[col], lens)})
    else:
        d.update({col: chainer(res[col])})                 
res2 = pd.DataFrame(d)
res2['GOterms'] = res2['GOterms'].str.strip() #remove whitespaces
#now add GOterm description and GOslims for each GOterm

def getGOdesc(GOterm):
    ind = GO.loc[GO['GO ID'] == GOterm].index
    if len(ind) > 0:
        ind = ind[0]
        GOdesc = GO.iloc[ind]['GO term']
        GOslim = GO.iloc[ind]['GO Slim(s)']
        return([GOdesc,GOslim])
    else:
        return(['',''])

res2['GOdesc'] = res2['GOterms'].map(lambda x:getGOdesc(x)[0])
res2['GOslims'] = res2['GOterms'].map(lambda x:getGOdesc(x)[1])

res2.to_csv(PT_output, header = True, sep = '\t', index = False)



    
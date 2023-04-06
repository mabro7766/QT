'''
script to generate regions file format required for the analyze_regions4.pl script
'''
#import modules

import pandas as pd
import sys 
import os
sys.path.append(os.path.abspath('Z:/mabro/QTs/paper/scripts/QT-GWAS/utilities')) #path to  where your config file is located containing all path variables (PT_<yourfile>)
from config_regions2 import *

data = pd.read_csv(PT_data, sep = '\t')

#add association number as the first column

data.insert(0, 'association number', range(1, len(data)+1))

#split the SNP column on '_' to separate the chromosome number from the location

data[['chromosome', 'location']] = data['SNP'].str.split('_', expand = True) 

#add a column containing the start and stop coordinates of the region (-10.000 to + 10.000 from the original location of the snp)
data['start'] = data['location'].astype('int64') - 10000
data['stop'] = data['location'].astype('int64') + 10000

#add a folder column
data['folder'] = 'folder'

#add trait number column

traits = data['Trait'].unique() #get all unique traits in a list
i = 0 # i will represent the trait number
d = {} #store trait and train number combinations in a dictionary
for trait in traits:
    i += 1
    d.update({trait:i})

data['traitnumber'] = data['Trait'].map(lambda x: d[x]) #obtain the matching trait number from the dictionary

#order columns and add 2 copies of the trait column to get to the regions format:
regions = data[['association number', 'chromosome', 'traitnumber', 'Trait', 'Trait', 'Trait', 'start', 'stop', 'folder']]
regions.to_csv(PT_output, sep = '\t', header = False, index = False) #write out the regions file as txt


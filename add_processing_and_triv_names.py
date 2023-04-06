# -*- coding: utf-8 -*-
"""
script to add processing and trivial names of the traits
Created on Wed Oct 20 16:39:39 2021

@author: mabro
"""
import pandas as pd
import os
import sys
sys.path.append(os.path.abspath('Z:/mabro/QTs/paper/scripts/QT-GWAS/utilities')) #path to  where your config file is located containing all path variables (PT_<yourfile>)
from config_add_processing_and_triv_names import *

#read in data
data = pd.read_csv(PT_data, sep = '\t')
procnames = pd.read_csv(PT_procnames, sep = '\t', header = None)
procnames.columns = ['name', 'procname']

#add processing names
def matchnames(trait, procnames):
    sel = procnames.loc[procnames['name'] == trait]['procname']#find the row in the procnames data with names and matching processing names of traits: the processing name
    #is then stored under the procname column
    if len(sel) == 0: #in case there are no matches return empty string
        procname = ''
    else:
        procname = sel.iloc[0] #else take the zeroth row to obtain the processing name
    return (procname)
data['processing name'] = data['trait'].map(lambda x: matchnames(x,procnames)) 

#add trivial names
chardata = pd.read_csv(PT_trivnames, sep = '\t')
def trivname(nodename):
    sel = chardata.loc[chardata['processing Name'] == nodename]
    if len(sel) > 0:
        return(sel.iloc[0]['trivial name'])
    else:
        return ('')
def compnr(nodename):
    sel = chardata.loc[chardata['processing Name'] == nodename]
    if len(sel) > 0:
        return(sel.iloc[0]['number'])
    else:
        return('')
data['trivial name'] = data['processing name'].map(lambda x: trivname(x))
data['compound number'] = data['processing name'].map(lambda x: compnr(x))
data.to_csv(PT_output, sep = '\t', index = False)

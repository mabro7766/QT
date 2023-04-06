"""
script for UMAP on SNP data, and comparison of subpopulations for QTs
a complete UMAP tutorial can be found here: https://umap-learn.readthedocs.io/en/latest/basic_usage.html
@author: mabro
"""

#import used modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import umap
from scipy import stats
import statsmodels.stats.multitest as multi
from matplotlib.colors import ListedColormap
from config_umap import * #config_umap.py in utilities folder contains all path to your file variables (PT_<yourfile>)

snpfile = pd.read_csv(PT_snpfile, sep = '\t') #read in SNP dataset with pandas
snpfile.insert(0, 'econr' ,  (range(1, len(snpfile) + 1))) #insert ecotype number
econrs = pd.read_csv(PT_ecodata, sep = '\t', header = None) #read in ecotype numbers with matching names, the ecotype numbers here match the ecotype numbers in the snp dataframe
econrs.columns = ['nr','name'] #rename the columns of the econrs dataframe
labels = list(econrs['name']) # get the labels (i.e. ecotype names) from the econrs dataframe

#UMAP on SNP data of ecotypes with longitude-latitude coordinates: ecotypes were sorted on longitude followed by latitude to
#simulate proximity: 8 ecotypes were ommited due to lack of coordinates (Kz-1, Da1-12,ICE138, ICE75, ICE112, ICE226, Tu-scha-9, Wal-Has-B4)

labels = pd.read_csv(PT_coords, sep = '\t', encoding = "ISO-8859-1") #read in ecotype coordinates data
econrs_todrop = labels[175:].nr #last 8 ecotypes of the df have no coordinates: retrieve their ecotype number to drop from the snp file
labels = labels[0:175] #drop last 8 rows with ecotypes without coordinates
labels['latitude'] = labels['latitude'].astype(float).round(2) #convert latitude to float with 2 decimals
labels['longitude'] = labels['longitude'].astype(float).round(2) #convert longitude to float with 2 decimals
labels['coord'] = labels['name'].astype(str)+ '-(' + labels['latitude'].astype(str) + ': ' + labels['longitude'].astype(str) + ')' #concatenate latitude and longiude into coordinates


#since we dropped the ecotypes lacking coordinates, we also have to drop them from snpfile on which we train the reducer
drp = []
for element in econrs_todrop.tolist():
    drp.append(element-1)
snpfile2 = snpfile.drop(drp)
reducer = umap.UMAP(random_state = 7766) #instantiate the class, with fixed random state for reproducibility 
# by fixin the random state be aware that this is a stochastic algorithm, you need to perform sufficient tests to 
# confirm that the main conclusions are not affected by this randomness
embedding = reducer.fit_transform(snpfile2) #train umap through fit, and immediately transform snp data to 2D dataset
embedding.shape
umap_df = pd.DataFrame(data = embedding, columns = ['dim1', 'dim2'], index = snpfile2['econr'].values) # index --> assign snpfile2 ecotype numbers as indices

#indices still in order of the original snp file: we want to match the new labels df order
umap_df['econr'] = snpfile2['econr'].values #used .values here otherwise the econrs get turned into floats

#sort umap_df according to the econr order in labels for plotting purpose
new_index = labels['nr'].tolist()
umap_sorted = umap_df.reindex(new_index)
names = labels['name'].tolist() 
umap_sorted['name'] = names #add the ecotype names to umap_sorted df


#plot UMAP results
labels_list = labels['coord'].tolist()
t = np.arange(175)
plt.figure(figsize=(12,14))
plt.scatter(umap_sorted['dim1'], umap_sorted['dim2'], c=t, cmap='viridis', s=40)
plt.tick_params(
    axis='both',          
    left = False,
    right = False,
    labelleft = False,
    bottom=False,      
    top=False,         
    labelbottom=False)
plt.gca().set_aspect('equal', 'datalim')
cbar = plt.colorbar(boundaries=np.arange(len(labels_list)+1)-0.5)
cbar.set_ticks(np.arange(len(labels_list)))
cbar.set_ticklabels(labels_list)
cbar.ax.tick_params(labelsize=6) 
plt.savefig(PT_output_umap,bbox_inches = 'tight', dpi = 300) #save the figure



#read in the metabolite data
abundances = pd.read_csv(PT_abundances,sep = '\t') #read in the abundances data, containing all metabolite abundances of all five replicates for each accession
QTs = pd.read_csv(PT_QTs,sep = '\t') #read in the QTs data, containing metabolites selected as QTs with p-value < 1x10E-06
new_header = abundances.iloc[0] #grab the first row for the header
abundances = abundances[1:] #take the data without the header row
abundances.columns = new_header #set the header row as the df header
abundances.rename(columns={np.nan:'metabolites'}, inplace=True) #rename the first column that contains the metabolites
QTs_list = QTs['QTs'].to_list()

#for each metabolite check in the five replicates of each ecotype if it is absent in more than half of the replicates
#currently the dataframe abundances_filtered ecotype columns are ordered alphabetical, but per repeat so we need to find the five repeats in the dataframe 
def checkabsence(metabolite): 
    producers = []
    nonproducers = []
    select = abundances.loc[abundances['metabolites'] == metabolite] #select the row matching the metabolite
    #cols = select[1:]
    #select = select[cols].apply(pd.to_numeric, errors='coerce', axis=1) #covert all columns except metabolite col to numeric values
    #for each ecotype we will retrieve the 5 columns from the abundance_filtered dataframe containing the five repeats
    for name in names: #names is the list of ecotypes that we defined above
        data = select.loc[:, select.columns.str[2:]==name]
        if int((data[data.columns.to_list()] == 0).sum(axis=1)) >= 3: #count how many repeats have a zero abundance, if >= 3 out of 5: then the accession is considered a nonproducer
            nonproducers.append(name) #add the ecotype to the nonproducers list
        else:
            producers.append(name) #otherwise the ecotype is producer: add to producers list
    return([producers, nonproducers]) #return list of producers and nonproducers

#plot UMAP colored according to producers vs non producers 
def plotumap(metabolite):    
    producers, nonproducers = checkabsence(metabolite) #call the producers and nonproducers for the specific metabolite
    eco_order = umap_sorted['name'].to_list()
    eco_labels = [] #here we will store 0 if ecotype is a non producer, and 1 if ecotype is a producer in the order of umap_sorted, for plotting
    for ecotype in eco_order:
        if ecotype in(producers):
            eco_labels.append(1)
        else:
            eco_labels.append(0)
    array_labels = np.array(eco_labels)
    colors = ListedColormap(['#a1dab4','#225ea8'])
    labels_list = ['nonproducer','producer']
    plt.figure(figsize=(12,14))
    plt.title("UMAP grouped producers and nonproducers of " + metabolite)
    scatter = plt.scatter(umap_sorted['dim1'], umap_sorted['dim2'], c=array_labels, cmap = colors, s=40)
    plt.tick_params(
        axis='both',          
        left = False,
        right = False,
        labelleft = False,
        bottom=False,      
        top=False,         
        labelbottom=False)
    plt.gca().set_aspect('equal', 'datalim')
    plt.legend(handles=scatter.legend_elements()[0], labels = labels_list)
    plt.savefig(PT_output_plots + "UMAP grouped producers and non producers of " + metabolite + ".png",bbox_inches = 'tight', dpi = 300) 
    plt.show()


#targeted for CYP450 associated metabolites
metabolite = 'M496.08_RT6.59'
plotumap(metabolite)
metabolite = 'M397.15_RT8.72'
plotumap(metabolite)
metabolite = 'M439.16_RT10.59'
plotumap(metabolite)
#N-acetyl isoleucin
metabolite = 'M172.24_RT5.89'
plotumap(metabolite)
#guanin benzoyl hexoside


#unknown associated to AT2G22960
metabolite = 'M463.15_RT2.48'
plotumap(metabolite)

#BGLU6 flavonoids
metabolite = 'M771.2_RT6.17'
plotumap(metabolite)
metabolite = 'M785.21_RT7.59'
plotumap(metabolite)
metabolite = 'M791.18_RT7.18'
plotumap(metabolite)

#GH3 jasmonate derivate
metabolite = 'M351.16_RT10.01'
plotumap(metabolite)

#guanin benzoyl hexoside
metabolite = 'M416.12_RT7.78'
plotumap(metabolite)
metabolite = 'M496.08_RT6.59'
plotumap(metabolite)

#Sulfated neolignans
metabolite = 'M469.08_RT6.54'
plotumap(metabolite)
metabolite = 'M499.09_RT9.62'
plotumap(metabolite)

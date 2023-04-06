'''
script to get the distributions of a trait in the subpopulations with/without snp

'''

def construct_nested_list(File):
    f = open(File,"r")
    content = f.read()
    f.close()
    content_list = content.split('\n') #split file per line 
    #del(content_list[0])  #delete header
    del(content_list[-1]) #delete last line (empty)
    #split each line per 'col' --> separated by '\t':
    nested_list = [] 
    for item in content_list:
        content_item = item.split('\t')
        nested_list.append(content_item)  #create nested list with each line as a list of the separate cols
    return(nested_list)

trtfile = construct_nested_list("//psb.ugent.be/research/groups/group_bioenergy/mabro/qualitative traits/abundances/ttrtf.txt") #contains data of traits: 0 or 1 for each ecotype (represented as 1-183): the cols are the ecotypes and the rows contain the trait values for each trait
snpfile = construct_nested_list("//psb.ugent.be/research/groups/group_bioenergy/mabro/qualitative traits/abundances/snps_sel.txt") #data of snps: 0 or 1 for each ecotype (represented as 1-183): the cols are the snps and the rows have the snp values for each ecotype --> transposed format of traitfile
ecotypes = construct_nested_list("//psb.ugent.be/research/groups/group_bioenergy/mabro/qualitative traits/abundances/ecotypenr_and_name.txt")
sign = construct_nested_list("//psb.ugent.be/research/groups/group_bioenergy/mabro/qualitative traits/abundances/significants.txt") #file with significant trait-snp pairs p < 10^-10
econr = construct_nested_list("//psb.ugent.be/research/groups/group_bioenergy/mabro/qualitative traits/abundances/econr.txt")
abundances = construct_nested_list("//psb.ugent.be/research/groups/group_bioenergy/mabro/qualitative traits/abundances/abundances.txt")
mgwasfile = construct_nested_list("//psb.ugent.be/research/groups/group_bioenergy/mabro/qualitative traits/abundances/mGWAST5_raw.txt")
significants = construct_nested_list("//psb.ugent.be/research/groups/group_bioenergy/mabro/qualitative traits/abundances/significants.txt")
final = construct_nested_list("//psb.ugent.be/research/groups/group_bioenergy/mabro/qualitative traits/abundances/final traits.txt")

#get the contingencytables                
def conttable(snpfile,snp,trtfile,trait):
    column = snpfile[0].index(snp) #find column at which the snp is
    searchlist_snp = [] 
    for line in snpfile[1:]:
        searchlist_snp.append(line[column]) #only use the col of the snp to count 0 and 1s
    import numpy as np
    arr = np.array(trtfile) #convert to numpy array
    row = list(arr[:,0]).index(trait) # slice out first col --> contains the traits and find the row at which current trait is
    searchlist_trait = list(arr[row,1:])# count in this row how many ones
    zero_zero = 0 # zero in snp and zero in trait val
    zero_one = 0 #zero in snp and 1 in trait
    one_zero = 0 #1 in snp and 0 in trait
    one_one = 0 #1 in both
    for i,snpval in enumerate(searchlist_snp):
        traitval = searchlist_trait[i]
        if snpval == "0" and traitval == "0": #ecotype does not have the snp nor the trait
            zero_zero += 1
        elif snpval == "0" and traitval == "1": #ecotype does not have the snp but does have the trait
            zero_one += 1
        elif snpval == "1" and traitval == "0": #ecotype has the snp but lacks the trait
            one_zero += 1
        elif snpval == "1" and traitval == "1": #ecotype has both snp and trait. 
            one_one += 1
    conttab = [[zero_zero,one_zero],[zero_one,one_one]]
    return(conttab) 


# first we need to devide the subpopulations: ecotypes with snp vs without
def subpop(snp):
    column = snpfile[0].index(snp) #find column at which the snp is
    searchlist_snp = []
    for line in snpfile[1:]:
        searchlist_snp.append(line[column]) #only use the col of the current snp
    snp_1 = []
    snp_0 = []
    for i,el in enumerate(searchlist_snp): #i = index in list --> i +1 is therefore the ecotype nr
        if el == "1":
            snp_1.append(str(i+1)) #snp_1 will be the subpop of ecotypes that have the snp
        elif el == "0":
            snp_0.append(str(i+1)) #snp_o will be the subpop of ecotypes that do not have the snp
    return ([snp_1,snp_0])

#for one trait: calculate average abundance across the 5 replicates of each ecotype    
def getaverage(abundances,trait):
    import numpy as np
    arr = np.array(abundances) #convert to numpy array
    row = list(arr[:,0]).index(trait) # slice out first col --> contains the traits and find the row at which current trait is
    traitlist = list(arr[row,1:])
    ecos = abundances[1][1:]
    d = {}
    for eco in ecos:
        if eco[2:5] != 'Col':
            indices = []
            for j,element in enumerate(ecos):
                if element[2:] == eco[2:]: # ecotype without replicate nr
                    indices.append(j)
            av = 0
            counter = 0
            for ind in indices:
                if traitlist[ind] != 'NA': 
                    counter += 1
                    av += float(traitlist[ind])
            av /= counter
            d.update({eco[2:]:av})
        else: #for col we have to make an exception because these accessions have random nrs after col
            indices = []
            for j,element in enumerate(ecos):
                if element[2:5] == 'Col':
                    indices.append(j)
            av = 0
            counter = 0
            for ind in indices:
                if traitlist[ind] != 'NA': 
                    counter += 1
                    av += float(traitlist[ind])
            av /= counter
            d.update({'Col':av})  
    return (d)
# sort dict by its values: sorted(d.items(), key=lambda x: x[1])

def getgwappnr(ecotype):
    l = construct_nested_list("//psb.ugent.be/shares/research/groups/group_bioenergy/mabro/qualitative traits/abundances/eco_gwappnrs.txt")
    for line in l:
        if line[0] == ecotype:
            return(line)
def getgwappname(nr):
    l = construct_nested_list("//psb.ugent.be/shares/research/groups/group_bioenergy/mabro/qualitative traits/abundances/eco_gwappnrs.txt")
    for line in l:
        if line[1] == nr:
            return(line)


'''
combine subpop and getaverage to get distributions per subpop
note in subpop the ecotypes are given by nrs, whereas in getaverage these are given by their names
so we will also need ecotypes variable in which the nrs are linked to the respective ecotype names
'''            
def av_subpop(snp,trait,abundances,ecotypes):
    subpop_1 = subpop(snp)[0] #this is the subpop for which the snp value is 1
    subpop_0 = subpop(snp)[1] #this is the subpop for which the snp value is 0
    av_subpop_1 = {}
    av_subpop_0 = {}
    d = getaverage(abundances,trait)
    for element in subpop_1:
        for el in ecotypes:
            if el[0] == element:
                ecotype = el[1] # remove the first two characters of el (=replicate nr)
                if ecotype[0:3] == 'Col': 
                    av = d['Col']
                    av_subpop_1.update({'Col':av})
                else:
                    av = d[ecotype]
                    av_subpop_1.update({ecotype:av})
    for element in subpop_0:
        for el in ecotypes:
            if el[0] == element:
                ecotype = el[1]# remove the first two characters of el (=replicate nr)
                if ecotype[0:3] == 'Col': 
                    av = d['Col']
                    av_subpop_0.update({'Col':av})
                else:
                    av = d[ecotype]
                    av_subpop_0.update({ecotype:av})
    return (av_subpop_1,av_subpop_0) # returns 2 dictionaries of the average abundance of a compound: for the population with snp 1 and snp 0 respectively



def getexceptions(d): #give in the subpopulation dictionary --> find which accessions are the exceptions
    storage = []
    for key in d.keys():
        zero = list(d.values()).count(0.0)
        if zero > len(list(d.values()))/2:# everything > 0 is an exception
            if d[key] != 0.0:
                storage.append([key,d[key]])
        else: # zeros are the exception
            if d[key] == 0.0:
                storage.append([key,d[key]])
    return(storage)

zeros = construct_nested_list('conttables_zero.txt')                
def newfishers(l): # get new fisher's exact tests for traits that have a zero in conttable --> replace by one
    import scipy.stats as stats
    for element in l:
        cont = []
        if element[6] == 0:
            cont.append(int('1'))
        else:
            cont.append(element[6])
            #also for 7,8,9
        

            
'''

def plot(d):
    import matplotlib.pyplot as plt
    plt.hist(d, bins=len(d))
    plt.show()'''

            
            

'''
check whether the associations can also be found in the mGWAS
'''
            
def assocmgwas(significants, mgwasfile):
    storage = []
    for line in significants:
        compound = line[0]
        snp = line[1]
        for line2 in mgwasfile:
            if line2[3] == compound:
                splitted = snp.split('_')
                if line2[1] == splitted[0] and int(line2[6]) <= int(splitted[1]) <= int(line2[7]):#chromosome nr matches and snp is located between 20 kb region
                    line.append('hit')
                    
        storage.append(line)
    return(storage)


'''
l2 = assocmgwas(significants, mgwasfile)
import csv
with open('mgwashits.csv','w',newline = '') as out: # write l2 as txt file
    csv_out=csv.writer(out)
    for row in l2:
        csv_out.writerow(row) '''

               
                
            
    
            
        
    
    

                
                    
            
                
                
                
        
                
        
    
    

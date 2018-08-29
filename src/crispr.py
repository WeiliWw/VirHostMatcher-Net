'''
This code generates feature values of CRISPR
'''

import os
import math
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
from .Variables import DB_HOST_CRISPR_PREFIX, TAXA_INFO
from numpy import isnan
'''
Preset variables and load tables
'''

db_host_crispr_prefix = DB_HOST_CRISPR_PREFIX 
# output_dir = INTERMEDIATE_RESULT 
taxa_info_file = TAXA_INFO 
#hash_table_file = HASH_TABLE 

db_host_crispr_prefix = os.path.expanduser(db_host_crispr_prefix)
#hash_table = pd.read_pickle(hash_table_file)

#dict_genome = hash_table.set_index(1).to_dict()[0]
taxa_info = pd.read_pickle(taxa_info_file)
taxa_info = taxa_info.set_index('hostNCBIName')

'''
Function: Run blast for a single query virus and return CRISPR signal
--------
Parameters:
    item: query file
    query_virus_dir
    numThreads
'''

def crisprSingle(item, query_virus_dir, output_dir, numThreads):
    query_name = item.split('.')[0]
    query_file = os.path.join(query_virus_dir, item)
    output_file = os.path.join(output_dir, query_name) + '.crispr'
    crispr_call = NcbiblastnCommandline(query=query_file,db=db_host_crispr_prefix,out=output_file,outfmt="'6 qacc sacc evalue'", evalue=1,gapopen=10,penalty=-1,
                                  gapextend=2,word_size=7,dust='no',
                                 task='blastn-short',perc_identity=90,num_threads=numThreads)
    crispr_call()
    '''
    Parse blast results
    '''
    if os.stat(output_file).st_size == 0:
        ind = False
        return ind, None
    else:
        query_res = pd.read_table(output_file,header = None)
        # Sanity check for blastn output format 
        query_res = query_res[query_res[1].apply(lambda x: x.count("|")) == 2]
        if query_res.shape[0] == 0:
            return False, None
        query_res[0] = query_name
        query_res[1] = query_res[1].apply(lambda x: x.split('|')[-2])
        #query_res[1] = [dict_genome[k] for k in list(query_res[1])]
        query_res[2] = -query_res[2].apply(math.log)
        df_crispr = query_res.groupby([0,1]).max().unstack(fill_value=0)
        ind = True
        return ind, df_crispr.set_index([[query_name]])

'''
Function: Gather available CRISPR signals
--------
Parameters:
    query_virus_dir
    numThreads
'''
def crispr_calculator(query_virus_dir, output_dir, numThreads):
    query_cont = []
    query_list = os.listdir(query_virus_dir)
    crispr_output_dir = os.path.join(output_dir, 'CRISPR/')
    try:
        os.stat(crispr_output_dir)
    except:
        os.mkdir(crispr_output_dir)
    for item in query_list:
        print('----Calculating crispr feature values for ',item,' ----')
        ind, df = crisprSingle(item, query_virus_dir, crispr_output_dir, numThreads)
        if ind:
            query_cont.append(df)
    print('----CRISPR intermediate files are stored in ',crispr_output_dir,' ----')
    if query_cont == []:
        return query_cont    # Return an empty list if no match for any queries
    else:
        df_concat = pd.concat(query_cont,axis =1,sort=False).groupby(axis=1,level=1,sort=False).sum().fillna(0)
        return df_concat
    
'''
Function: Augment signals to the entire matrix and generate values at genus level
--------
Parameters:
    df_input
    virus_index: query index
    host_index
'''

def uniGenus(df_input, virus_index, host_index):
    df_pseudo = pd.DataFrame(index=virus_index, columns=host_index).fillna(0)
    df_full = pd.concat([df_pseudo, df_input],sort=False).groupby(level=0,sort=False).sum().fillna(0)
    df_full = df_full.loc[virus_index][host_index]
    df_full.loc['hostGenus'] = taxa_info.loc[host_index]['hostGenus']
    dict_genera = {}
    for i in df_full:
        col = df_full[i].rename(None)
        genus = col['hostGenus']
        if type(genus) is not str:       # not a genus name
            if isnan(genus): continue
        if genus in dict_genera:
            dict_genera[genus] = pd.concat([dict_genera[genus], col],axis=1,sort=False).max(axis=1)
        else: dict_genera[genus] = col
    for i in df_full:
        genus = df_full[i]['hostGenus']
        if type(genus) is not str:
            if isnan(genus): continue
        df_full[i] = dict_genera[genus] 
    return df_full.loc[virus_index,:]
'''
def uniGenus(df_input, host_index):
    virus_index = df_input.index
    df_pseudo = pd.DataFrame(index=virus_index, columns=host_index).fillna(0)
    df_full = pd.concat([df_pseudo, df_input]).groupby(level=0).sum().fillna(0)
    df_full = df_full.loc[virus_index][host_index]
    df_full.loc['hostGenus'] = taxa_info.loc[host_index]['hostGenus']
    for genera in taxa_info.loc[host_index]['hostGenus'].unique():
        if pd.notnull(genera):
            idx = (df_full.loc['hostGenus'] == genera)
            df_full.loc[virus_index,idx] = df_full.loc[virus_index,idx].max(axis=1)
    return df_full.loc[virus_index,:]
'''








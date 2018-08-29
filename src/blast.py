'''
This code generates feature values of BLAST (fraction of aligned virus 
genome)
'''

import os
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
from .Variables import DB_HOST_PREFIX, HASH_TABLE

'''
Preset variables and load hash table data
'''
db_host_prefix = DB_HOST_PREFIX 
# output_dir = INTERMEDIATE_RESULT 
hash_table_file = HASH_TABLE 

db_host_prefix = os.path.expanduser(db_host_prefix)
hash_table = pd.read_pickle(hash_table_file)
hash_table[1] = hash_table[1].apply(lambda x: x.split('.')[0])
dict_genome = hash_table.set_index(1).to_dict()[0]


'''
Auxiliary function
'''
'''
def cal_perc(x):
    indicator = [0]*x[4]
    for i in range(len(x[2])):
        indicator[x[2][i]:x[3][i]] = [1]*(x[3][i]-x[2][i]+1)
    return sum(indicator)/x[4]
'''
def cal_perc(x):
    indicator = [0]*x[4]
    for i in range(len(x[2])):
        indicator[x[2][i]:x[3][i]] = [1]*(x[3][i]-x[2][i]+1)
    return sum(indicator)
'''
Function: Run blast for a single query file and return the mapping percentage
--------
    item: query file
    query_virus_dir
    numThreads
'''
'''
def blastSingle(item, query_virus_dir, numThreads):
    query_name = item.split('.fa')[0]
    query_file = os.path.join(query_virus_dir, item)
    output_file = os.path.join(output_dir, query_name) + '.blast'
    blast_call = NcbiblastnCommandline(query=query_file,db=db_host_prefix, 
                                       out=output_file,outfmt="'6 qacc sacc qstart qend qlen'", 
                                       evalue=0.01,gapopen=10,penalty=-1,reward=1,gapextend=2,
                                       word_size=11,perc_identity=90,num_threads=numThreads)
    blast_call()
    if os.stat(output_file).st_size == 0:
        ind = False
        return ind, None
    else:
        query_res = pd.read_table(output_file,header = None)
        # need to make sure a same value for the last column
        # map headers to genome names
        #query_res[0] = query_name
        query_res[1] = [dict_genome[k] for k in list(query_res[1])]
        df_blast_positions = query_res.groupby([0,1]).agg({2: lambda x: tuple(x-1), 
                                                           3: lambda x: tuple(x-1), 4: min})
        df_blast = df_blast_positions.apply(lambda x: cal_perc(x), 
                                            axis = 1).unstack(fill_value=0)
        ind = True
        return ind, df_blast.set_index([[query_name]])
'''
def blastSingle(item, query_virus_dir, output_dir, numThreads):
    query_name = item.split('.')[0]
    query_file = os.path.join(query_virus_dir, item)
    output_file = os.path.join(output_dir, query_name) + '.blast'
    blast_call = NcbiblastnCommandline(query=query_file,db=db_host_prefix, 
                                       out=output_file,outfmt="'6 qacc sacc qstart qend qlen'", 
                                       evalue=0.01,gapopen=10,penalty=-1,reward=1,gapextend=2,
                                       word_size=11,perc_identity=90,num_threads=numThreads)
    blast_call()

    '''
    Parse blast results for a single file
    '''
    if os.stat(output_file).st_size == 0:
        ind = False
        return ind, None
    else:
        with open(query_file) as f:
            query_len = len(f.read())  # bp
        query_res = pd.read_table(output_file,header = None)
        # need to make sure a same value for the last column
        # map headers to genome names
        #query_res[0] = query_name
        query_res[1] = [dict_genome[k] for k in list(query_res[1])]
        df_blast_positions = query_res.groupby([0,1]).agg({2: lambda x: tuple(x-1), 
                                                           3: lambda x: tuple(x-1), 4: min})
        df_blast_positions.index = df_blast_positions.index.droplevel()
        df_blast_perc = df_blast_positions.apply(lambda x: cal_perc(x), 
                                                 axis = 1)/query_len
        sr_blast = df_blast_perc.groupby(level=0,sort=False).apply(sum)
        ind = True
        return ind, pd.DataFrame({query_name: sr_blast}).T
    
    
    
'''
Calculate blast feature values for all query viruses
Return a df with all available maps (not a FULL matrix of ALL hosts and queries)
--------
Parameters:
    query_virus_dir
    virus_index
    host_index
    numThreads
'''
    
def blast_calculator(query_virus_dir, virus_index, host_index, output_dir, numThreads=1):
    blast_output_dir = os.path.join(output_dir, 'BLAST/')
    try:
        os.stat(blast_output_dir)
    except:
        os.mkdir(blast_output_dir)
    query_cont = []
    query_list = os.listdir(query_virus_dir)
    for item in query_list:
        print('----Calculating blast feature values for ',item,' ----')
        ind, df = blastSingle(item, query_virus_dir, blast_output_dir, numThreads)
        if ind:
            query_cont.append(df)
    df_pseudo = pd.DataFrame(index=virus_index, columns=host_index).fillna(0)
    if query_cont == []:
        print('----BLAST intermediate files are stored in ',blast_output_dir,' ----')
        return df_pseudo   # return a zero matrix if no match
    else:
        tmp = pd.concat(query_cont, axis =1, sort=False).groupby(axis=1,level=0,sort=False).sum()
        print('----BLAST intermediate files are stored in ',blast_output_dir,' ----')
        return pd.concat([tmp, df_pseudo],sort=False).groupby(level=0,sort=False).sum().fillna(0).loc[virus_index][host_index] 
    
    
    
 

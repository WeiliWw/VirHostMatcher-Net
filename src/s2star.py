'''
Generates necessary results of s2star 
'''


import os
import pandas as pd
from . import tools
from .Variables import *
'''
Preset variables
### to be checked!
'''

original_interaction_file = S2STAR_BENCH_HOST #'./prepare_data/relationMat_352by31986.csv'
hash_dir = D2STAR_HASH #'/home/rcf-40/weiliw/panasas/v-h-SAGs/d2star_dissimilarity/hash'
output_dir = INTERMEDIATE_RESULT #'./intermediate_res/'
pseudo_host = PSEUDO_HOST #'./prepare_data/pseudo_host/'
pseudo_virus = PSEUDO_VIRUS #'./prepare_data/pseudo_virus/'

'''
Function: Calculate distances within query viruses
--------
Parameters:
    query_virus_dir = './test_query/'

'''
def intra_query(query_virus_dir):
    print('----Start calculating s2* part I. This may take a few seconds...----')
    tools.cafe(['CAFE','-F1',query_virus_dir, '-F2', query_virus_dir, 
                '-K', '6', '-M', '2', '-D', 'D2star', '-S', hash_dir, 
                '-R', '-O', output_dir + 'intra_query_virus'])
    print('----Finished calculating s2* part I----')
    # Read results and return the df
    results_file = output_dir + 'intra_query_virus.D2star.plain'
    results = pd.read_csv(results_file, index_col=0)
    return 1 - 2*results

'''
Function: Calculate distances between query viruses and 352 'benchmark' 
viruses
--------
Parameters:
    query_virus_dir

'''
def query_bench(query_virus_dir):
    print('----Start calculating s2* part II. This may take a few minutes...----')
    tools.cafe(['CAFE','-F1',query_virus_dir, '-F2', pseudo_virus, 
                '-K', '6', '-M', '2', '-D', 'D2star', '-S', hash_dir, 
                '-R', '-O', output_dir + 'between_query_virus_bench'])
    print('----Finished calculating s2* part II----')
    # Read results and return the df
    results_file = output_dir + 'between_query_virus_bench.D2star.plain'
    results = pd.read_csv(results_file, index_col=0)
    return 1 - 2*results


'''
Function: Calculate distances between query viruses and candidate hosts
--------
Parameters:
    query_virus_dir

'''
def query_host(query_virus_dir):
    print('----Start calculating s2* part III. This may take up to a few hours...----')
    tools.cafe(['CAFE','-F1',query_virus_dir, '-F2', pseudo_host, 
                '-K', '6', '-M', '2', '-D', 'D2star', '-S', hash_dir, 
                '-R', '-O', output_dir + 'inter_query_virus_hosts'])
    print('----Finished calculating s2* part III----')
    # Read results and return the df
    results_file = output_dir + 'inter_query_virus_hosts.D2star.plain'
    results = pd.read_csv(results_file, index_col=0)
    return 1 - 2*results


'''
Function: Combine previous d2star results and return virus-virus and virus-host
distance matrices and interaction matrix
--------
Parameters:
    query_virus_dir

'''

'''
def s2star_caclculator(query_virus_dir):
    #mat_original_v_h = pd.read_csv(original_s2star_virus_host_file, index_col=0)
    #mat_intra_bench = pd.read_csv(original_bench_file, index_col=0)
    mat_original_interaction = pd.read_csv(original_interaction_file, index_col=0)
    print('----Calculation of s2* is split into three parts----')
    mat_intra_query = intra_query(query_virus_dir)
    mat_query_bench = query_bench(query_virus_dir)
    s2star_query_host = query_host(query_virus_dir)
    ################ Add an assertion ? ############
    #s2star_virus_host = pd.concat([mat_original_v_h, mat_query_host])
    #s2star_intra_virus = pd.concat([mat_intra_bench, mat_query_bench, 
    #                               mat_query_bench.T,
    #                               mat_intra_query]).groupby(level=0).sum()
    s2star_query_virus = pd.concat([mat_query_bench,
                                    mat_intra_query]).groupby(level=0).sum()
    host_index = s2star_query_host.columns
    virus_index = s2star_query_virus.columns
    query_index = mat_intra_query.index
    # rearrange df
    s2star_query_virus = s2star_query_virus.loc[query_index]
    s2star_query_host = s2star_query_host.loc[query_index]
    pseudo_interaction = pd.DataFrame(index=virus_index, columns=host_index).fillna(0)
    df_interaction = pd.concat([mat_original_interaction,pseudo_interaction]).groupby(level=0).sum().fillna(0)
    df_interaction = df_interaction.loc[virus_index][host_index]
    return s2star_query_host, s2star_query_virus, df_interaction
'''

def s2star_caclculator(query_virus_dir, ifShort):
    #mat_original_v_h = pd.read_csv(original_s2star_virus_host_file, index_col=0)
    #mat_intra_bench = pd.read_csv(original_bench_file, index_col=0)
    if ifShort:
        mat_original_interaction = pd.read_csv(original_interaction_file, index_col=0)
        print('----Calculation of s2* is split into two parts----')
        mat_intra_query = intra_query(query_virus_dir)
        mat_query_bench = query_bench(query_virus_dir)
        host_index = mat_original_interaction.columns
        query_index = mat_intra_query.index
        virus_index = mat_original_interaction.index.union(query_index)
        s2star_query_virus = pd.concat([mat_query_bench,
                                        mat_intra_query]).groupby(level=0).sum()
        # rearrange df
        s2star_query_host = []
        s2star_query_virus = s2star_query_virus.loc[query_index, virus_index]
        pseudo_interaction = pd.DataFrame(index=virus_index, columns=host_index).fillna(0)
        df_interaction = pd.concat([mat_original_interaction,pseudo_interaction]).groupby(level=0).sum().fillna(0)
        df_interaction = df_interaction.loc[virus_index][host_index]
    else:
        mat_original_interaction = pd.read_csv(original_interaction_file, index_col=0)
        print('----Calculation of s2* is split into three parts----')
        mat_intra_query = intra_query(query_virus_dir)
        mat_query_bench = query_bench(query_virus_dir)
        host_index = mat_original_interaction.columns
        query_index = mat_intra_query.index
        virus_index = mat_original_interaction.index.union(query_index)
        s2star_query_virus = pd.concat([mat_query_bench,
                                        mat_intra_query]).groupby(level=0).sum()
        # s2star feature
        s2star_query_host = query_host(query_virus_dir)    
        # rearrange df
        s2star_query_host = s2star_query_host.loc[query_index, host_index]
        s2star_query_virus = s2star_query_virus.loc[query_index, virus_index]
        pseudo_interaction = pd.DataFrame(index=virus_index, columns=host_index).fillna(0)
        df_interaction = pd.concat([mat_original_interaction,pseudo_interaction]).groupby(level=0).sum().fillna(0)
        df_interaction = df_interaction.loc[virus_index][host_index]
    return s2star_query_host, s2star_query_virus, df_interaction
















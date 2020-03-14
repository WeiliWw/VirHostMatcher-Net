'''
Generates necessary results of s2star 
'''
import os, sys, pickle
import pandas as pd
#from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
#from . import tools
from .tools import kmer_count
from .Variables import TABLE_HOST, TABLE_BENCH, TABLE_INTER
'''
Preset variables
### to be checked!
'''
# output_dir = INTERMEDIATE_RESULT #'./intermediate_res/'
# tables = TABLES
# original_interaction_file = S2STAR_BENCH_HOST 
# hash_dir = D2STAR_HASH #'/home/rcf-40/weiliw/panasas/v-h-
# pseudo_host = PSEUDO_HOST #'./prepare_data/pseudo_host/'
# pseudo_virus = PSEUDO_VIRUS #'./prepare_data/pseudo_virus/'

'''
Helper functions for calculating s2star
'''
Alphabeta = ['A', 'C', 'G', 'T']
Alpha_dict = dict(zip(Alphabeta, range(4)))

def get_transition(count_array):
    shape = len(count_array)
    transition_array = count_array.reshape(shape//4, 4)
    with np.errstate(divide='ignore', invalid='ignore'):
        transition_array = (transition_array / np.sum(transition_array, 1)[:, np.newaxis])
        transition_array[np.isnan(transition_array)] = 0
    return transition_array


def get_expect(a_M_count, a_trans, K, M):
    a_expect = a_M_count
    for _ in range(K-M):
        a_trans = np.tile(a_trans, (4, 1))
        a_expect = (a_expect[:,np.newaxis] * a_trans).flatten()
    return a_expect

def get_f(a_K_count, a_expect):
    with np.errstate(divide='ignore', invalid='ignore'):
        f = (a_K_count-a_expect)/np.sqrt(a_expect)
        f[np.isnan(f)]=0
    return f

def get_all_f(Dir, K, order, Reverse, numThreads):
    M = order + 1
    sequence_list = os.listdir(Dir)
    f_matrix = np.ones((len(sequence_list), 4**K))
    flag = False
    for i, seq in enumerate(sequence_list):
        seqfile = os.path.join(Dir, seq)
        K_count = np.array(kmer_count(seqfile, numThreads, Reverse, K))
        if K_count[0] == -1:            # sanity check for the fasta files
            print('The query file {} contains invalid chars, please make sure it is a valid fasta file.'.format(seq))
            flag = True
        if np.sum(K_count) == 0:
            print('The query file {} is empty, please double check the file.'.format(seq))
            flag = True
        M_count = np.array(kmer_count(seqfile, numThreads, Reverse, M))
        trans = get_transition(M_count)
        expect = get_expect(M_count, trans, K, M)
        f_matrix[i] = get_f(K_count, expect)
    name_list = [x.rsplit('.', 1)[0] for x in sequence_list]
    if flag: sys.exit('Program terminated. Please check error info above.')
    return f_matrix, name_list

def cosine_similarity(f1, f2):
    n1 = np.linalg.norm(f1,axis=1,keepdims=True)
    n2 = np.linalg.norm(f2,axis=1,keepdims=True)
    norms = np.dot(f1,f2.T)
    prod = np.dot(f1,f2.T)
    return prod/(np.dot(n1,n2.T))


'''
Function: Calculate distances within query viruses
--------
Parameters:
    query_virus_dir = './test_query/'

'''



def s2_query(query_virus_dir, ifShort, numThreads):
    print('----Start calculating s2* part I... ----')
    virus_mat, virus_list = get_all_f(query_virus_dir, 6, 2, True, numThreads)
    # s2_intra_mat = cosine_similarity(virus_mat, virus_mat)
    print('----Finished calculating s2* part I----')
    print('----Start calculating s2* part II... ----')
    ## read bench file
    bench_mat, bench_list = pickle.load(open(TABLE_BENCH, 'rb'))
    s2_virus_bench_mat = cosine_similarity(virus_mat, bench_mat)
    print('----Finished calculating s2* part II----')
    if ifShort:
        return pd.DataFrame(s2_virus_bench_mat, index=virus_list, columns=bench_list), None
    else:
        ## read host mat file
        print('----Start calculating s2* part III... ----')
        host_mat, host_list = pickle.load(open(TABLE_HOST, 'rb'))
        s2_host_virus_mat = cosine_similarity(virus_mat, host_mat)
        print('----Finished calculating s2* part III----')
        ## return three matrices
        return pd.DataFrame(s2_virus_bench_mat, index=virus_list, columns=bench_list), pd.DataFrame(s2_host_virus_mat, index=virus_list, columns=host_list)



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
'''
Function: Calculate distances between query viruses and 352 'benchmark' 
viruses
--------
Parameters:
    query_virus_dir

'''


# def query_bench(query_virus_dir):
#     print('----Start calculating s2* part II. This may take a few minutes...----')
#     tools.cafe(['CAFE','-F1',query_virus_dir, '-F2', pseudo_virus, 
#                 '-K', '6', '-M', '2', '-D', 'D2star', '-S', hash_dir, 
#                 '-R', '-O', output_dir + 'between_query_virus_bench'])
#     print('----Finished calculating s2* part II----')
#     # Read results and return the df
#     results_file = output_dir + 'between_query_virus_bench.D2star.plain'
#     results = pd.read_csv(results_file, index_col=0)
#     return 1 - 2*results


'''
Function: Calculate distances between query viruses and candidate hosts
--------
Parameters:
    query_virus_dir

'''
# def query_host(query_virus_dir):
#     print('----Start calculating s2* part III. This may take up to a few hours...----')
#     tools.cafe(['CAFE','-F1',query_virus_dir, '-F2', pseudo_host, 
#                 '-K', '6', '-M', '2', '-D', 'D2star', '-S', hash_dir, 
#                 '-R', '-O', output_dir + 'inter_query_virus_hosts'])
#     print('----Finished calculating s2* part III----')
#     # Read results and return the df
#     results_file = output_dir + 'inter_query_virus_hosts.D2star.plain'
#     results = pd.read_csv(results_file, index_col=0)
#     return 1 - 2*results


'''
Function: Combine previous d2star results and return virus-virus and virus-host
distance matrices and interaction matrix
--------
Parameters:
    query_virus_dir

'''

#
#def s2star_caclculator(query_virus_dir, ifShort, numThreads):
#    mat_original_interaction = pd.read_hdf(tables,'interaction') # 352 by 31k
#    if ifShort:
#        print('----Calculation of s2* is split into two parts----')
#    else:
#        print('----Calculation of s2* is split into three parts----')
#    mat_intra_query, mat_query_bench, s2star_query_host = s2_query(query_virus_dir, ifShort, numThreads)
##         mat_intra_query = intra_query(query_virus_dir)
##         mat_query_bench = query_bench(query_virus_dir)
#    host_index = mat_original_interaction.columns
#    query_index = mat_intra_query.index
#    virus_index = mat_original_interaction.index.union(query_index)
#    s2star_query_virus = pd.concat([mat_query_bench,
#                                    mat_intra_query],sort=False).groupby(level=0,sort=False).sum()
#    # rearrange df
#    if not ifShort:
#        # rearrange df
#        s2star_query_host = s2star_query_host.loc[query_index, host_index]
#    s2star_query_virus = s2star_query_virus.loc[query_index, virus_index]
#    pseudo_interaction = pd.DataFrame(index=virus_index, columns=host_index).fillna(0)
#    df_interaction = pd.concat([mat_original_interaction,pseudo_interaction]).groupby(level=0).sum().fillna(0)
#    df_interaction = df_interaction.loc[virus_index][host_index]
#    return s2star_query_host, s2star_query_virus, df_interaction
#
def s2star_caclculator(query_virus_dir, ifShort, numThreads):
    mat_original_interaction = pd.read_csv(TABLE_INTER) # 352 by 31k
    if ifShort:
        print('----Calculation of s2* is split into two parts----')
    else:
        print('----Calculation of s2* is split into three parts----')
    mat_query_bench, s2star_query_host = s2_query(query_virus_dir, ifShort, numThreads)
    host_index = mat_original_interaction.columns
    query_index = mat_query_bench.index
    virus_index = query_index.tolist() + mat_original_interaction.index.tolist()
    #s2star_query_virus = pd.concat([mat_query_bench,
    #                                mat_intra_query], sort=False).groupby(level=0, sort=False).sum().loc[query_index, virus_index]
    if not ifShort:
        s2star_query_host = s2star_query_host.loc[query_index, host_index]
    #s2star_query_virus = s2star_query_virus.loc[query_index, virus_index]
    #pseudo_interaction = pd.DataFrame(index=virus_index, columns=host_index).fillna(0)
    #df_interaction = pd.concat([mat_original_interaction,pseudo_interaction]).groupby(level=0).sum().fillna(0)
    #df_interaction = df_interaction.loc[virus_index][host_index]
    #s2star_query_virus.values[[np.arange(s2star_query_virus.shape[0])]*2] = 0
    return s2star_query_host, mat_query_bench, mat_original_interaction 



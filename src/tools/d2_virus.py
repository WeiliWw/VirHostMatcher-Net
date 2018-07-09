
from _count import kmer_count
from scipy.spatial.distance import cosine
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import time
import os


Alphabeta = ['A', 'C', 'G', 'T']
Alpha_dict = dict(zip(Alphabeta, range(4)))
# 
# def num2nuc(num, k):
#     num_bin = format(num, '0%sb'%(k*2))
#     return ''.join([Alphabeta[int(num_bin[i:i+2], 2)] for i in range(0, len(num_bin), 2)])
# 
# def shift(num, K, M, x):
#     mask = 2**(2*M)-1
#     num = num >> ((K-M)*2 - 2*x)
#     num &= mask
#     return num
# 
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
   
def get_all_f(Dir, K, M, Reverse, Num_Threads):
    sequence_list = os.listdir(Dir)[:1000]
    f_matrix = np.ones((len(sequence_list), 4**K))
    for i, seq in enumerate(sequence_list):
        seqfile = os.path.join(Dir, seq)
        K_count = np.array(kmer_count(seqfile, K, Reverse, Num_Threads))
        M_count = np.array(kmer_count(seqfile, M, Reverse, Num_Threads))
        trans = get_transition(M_count)
        expect = get_expect(M_count, trans, K, M)
        f_matrix[i] = get_f(K_count, expect)
    return f_matrix
        
    
if __name__ == "__main__":
    begin = time.time()
    Dir1 = "/home/weili/Documents/datasets/EVGs"
    Dir2 = "/home/weili/Documents/datasets/EVGs"
    Num_Threads = 2 
    K = 6
    M = 3
    Reverse = True
    f1 = get_all_f(Dir1, K, M, Reverse, Num_Threads)
    f2 = get_all_f(Dir2, K, M, Reverse, Num_Threads)
    d2_matrix = 0.5 * (1 - cosine_similarity(f1, f2))
    np.fill_diagonal(d2_matrix, 0)
    print(time.time()-begin)

'''
This code generates likelihood generated from WIsH 
'''
import os
import pandas as pd
from . import tools 
from .Variables import WISH_HOST_MODELS

'''
Preset variables 
'''
host_model_dir = WISH_HOST_MODELS 
# output_dir = INTERMEDIATE_RESULT 

'''
Function: Fit query viruses into host models and return likelihood matrix
--------
Parameters:
    query_virus_dir
    virus_index
    host_index
'''

def wish_llkd_calculator(query_virus_dir, virus_index, host_index, output_dir, numThreads):
#def wish_llkd_calculator(query_virus_dir, numThreads, output_dir,  host_model_dir = WISH_HOST_MODELS):
# import src.wish
# src.wish.wish_llkd_calculator('test_query/', 8, 'tmp', '/home/rcf-40/weiliw/panasas/v-h-WIsH/host2695_model/')
    print('----Fitting models in WIsH...----')
    tools.wish(query_virus_dir, host_model_dir, output_dir, 'predict', numThreads)
    print('----WIsH calculation finished.----')
    llkh = pd.read_table(os.path.join(output_dir, 'llikelihood.matrix'))
    return llkh.T.loc[virus_index][host_index]
'''
def wish_llkd_calculator(query_virus_dir, virus_index, host_index, numThreads):
    print('----Fitting models in WIsH. This may take up to a few hours...----')
    tools.wish(['wish','-c','predict','-m',host_model_dir,'-g',query_virus_dir,
                '-r',output_dir,'-b', '-t',str(numThreads)])
    print(' WIsH finished. ')
    llkh = pd.read_table(os.path.join(output_dir, 'llikelihood.matrix'))
    return llkh.T.loc[virus_index][host_index]
'''

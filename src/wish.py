'''
This code generates feature values of WIsH likelihood 
'''
import os
import pandas as pd
from . import tools 
from .Variables import *

'''
Preset variables 
'''
host_model_dir = WISH_HOST_MODELS 
output_dir = INTERMEDIATE_RESULT 

'''
Function: Fit query viruses into host models and return likelihood matrix
--------
Parameters:
    query_virus_dir
    virus_index
    host_index
'''
def wish_llkd_calculator(query_virus_dir, virus_index, host_index):
    print('----Fitting models in WIsH. This may take up to a few hours...----')
    tools.wish(['wish','-c','predict','-g',query_virus_dir,'-m',host_model_dir,
                '-r',output_dir,'-t','12','-b'])
    llkh = pd.read_table(os.path.join(output_dir, 'llikelihood.matrix'))
    return llkh.T.loc[virus_index][host_index]


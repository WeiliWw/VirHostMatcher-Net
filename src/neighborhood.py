'''
This code is used to generate positive and negative SV's, given interaction matrix,
intra virus and virus-host s2star matrices
'''

import pandas as pd

def neighborhood_calculator(df_query_virus, df_interaction):
    '''
    ####df_query_host: s2star matrix: (Query Viruses) * Hosts
    df_query_virus: s2star matrix: Q * (B + Q)
    df_interaction: binary matrix: (Bench + Query Viruses) * Hosts
    '''
    print('----Start calculating network neighborhood feature values...----')
    pos_interaction = df_interaction.apply(lambda x: x/sum(x) ,axis=0).fillna(0)
    pos_SV = df_query_virus.dot(pos_interaction)
    neg_interaction = (1 - df_interaction).apply(lambda x: x/sum(x),
                                                 axis=0).fillna(0)
    neg_SV = df_query_virus.dot(neg_interaction)
    print('----Finished Calculating network neighborhood feature values----')
    return pos_SV, neg_SV


'''
Build a class for queries
'''
import src.s2star
import src.crispr
import src.neighborhood
import src.wish
import src.blast
from src.Variables import REGRESSION_COEFFICIENTS, REGRESSION_COEFFICIENTS_SHORT, TAXA_INFO
import pandas as pd
import numpy as np

class HostPredictor:
    def __init__(self, query_virus_dir, ifShort=False, numThreads=1):
        '''
        Attributes calculated:
            s2star, posSV, negSV, crispr, blast, wish
        '''
        self._short = ifShort
        self.s2star, self._query_virus, self._df_interaction = src.s2star.s2star_caclculator(query_virus_dir, ifShort)
        self.posSV, self.negSV = src.neighborhood.neighborhood_calculator(self._query_virus, self._df_interaction)
        self._virus_index = self._query_virus.index  # query virus index
        self._host_index = self._df_interaction.columns
        self._crispr_signals = src.crispr.crispr_calculator(query_virus_dir, numThreads)
        self.crispr = src.crispr.uniGenus(self._crispr_signals, self._virus_index, self._host_index)
        # wish , if statement?
        self.blast = src.blast.blast_calculator(query_virus_dir, self._virus_index, self._host_index, numThreads)
        if ifShort:
            self.wish = src.wish.wish_llkd_calculator(query_virus_dir, self._virus_index, self._host_index)
        else:
            self.wish = None
        #self.blast = src.blast.blast_calculator(query_virus_dir, self._virus_index, self._host_index, numThreads)
            
    def getScores(self):
        '''
        Output a prediction table
        '''
        #pred = pd.DataFrame(index=virus_index, columns=host_index).fill
        if self._short:
            score = REGRESSION_COEFFICIENTS_SHORT[0]*pd.DataFrame(index=self._virus_index, columns=self._host_index).fillna(1) +  \
            REGRESSION_COEFFICIENTS_SHORT[1]*self.wish +  \
            REGRESSION_COEFFICIENTS_SHORT[2]*self.posSV +  \
            REGRESSION_COEFFICIENTS_SHORT[3]*self.negSV +  \
            REGRESSION_COEFFICIENTS_SHORT[4]*self.crispr +  \
            REGRESSION_COEFFICIENTS_SHORT[5]*self.blast
        else:
            score = REGRESSION_COEFFICIENTS_SHORT[0]*pd.DataFrame(index=self._virus_index, columns=self._host_index).fillna(1) +  \
            REGRESSION_COEFFICIENTS_SHORT[1]*self.s2star +  \
            REGRESSION_COEFFICIENTS_SHORT[2]*self.posSV +  \
            REGRESSION_COEFFICIENTS_SHORT[3]*self.negSV +  \
            REGRESSION_COEFFICIENTS_SHORT[4]*self.crispr +  \
            REGRESSION_COEFFICIENTS_SHORT[5]*self.blast
        self.score = 1- 1/(pd.DataFrame(score, dtype=np.float).apply(np.exp)+1)
        
    def prediction(self, topN):   # return a dictionary of queries
        dict_pred = {}
        # read in taxa info:
        taxa_info = pd.read_pickle(TAXA_INFO)
        taxa_info = taxa_info.set_index('hostNCBIName')
        '''
        If the highest scores have a tie more than topN, output them all
        '''
        for query in self.score.index:
            query_score = self.score.loc[query]
            topIdx = list()
            topNum = 0
            while topNum < topN and topNum <= len(self._host_index):
                idx = query_score[query_score == query_score.max()].index
                topIdx.extend(idx)
                topNum = len(topIdx)
                query_score = query_score.drop(idx)
                #print("size of query: ", query_score.size)
            taxa = taxa_info.loc[topIdx]
            taxa['score'] = self.score.loc[query][topIdx]
            dict_pred[query] = taxa
        return dict_pred
            
            
        
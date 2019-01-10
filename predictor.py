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
import os

class HostPredictor:
    def __init__(self, query_virus_dir, ifShort, intermediate_dir, genome_list, numThreads=1):
        '''
        Attributes calculated:
            s2star, posSV, negSV, crispr, blast, wish
        '''
        try:
            os.stat(intermediate_dir)
        except:
            os.mkdir(intermediate_dir)
        print("Intermediate results will be stored in ", intermediate_dir)
        
        self._short = ifShort
        self.s2star, self._query_virus, self._df_interaction = src.s2star.s2star_caclculator(query_virus_dir, ifShort, numThreads)
        self.posSV, self.negSV = src.neighborhood.neighborhood_calculator(self._query_virus, self._df_interaction)
        self._virus_index = self._query_virus.index  # query virus index
        if genome_list is None:                     # when no genome list specified
            self._host_index = self._df_interaction.columns
            genomes = None
        else:
            with open(genome_list) as f:
                genomes = f.readlines()
            genomes = [i.rstrip() for i in genomes]
            self.posSV = self.posSV[genomes]
            self.negSV = self.negSV[genomes]
            if not ifShort:
                self.s2star = self.s2star[genomes]
            self._host_index = genomes
        if ifShort:
            self.wish = src.wish.wish_llkd_calculator(query_virus_dir, self._virus_index, self._host_index, intermediate_dir, numThreads)
        else:
            self.wish = None
        self.blast = src.blast.blast_calculator(query_virus_dir, self._virus_index, self._host_index, genomes, intermediate_dir, numThreads)
        self._crispr_signals = src.crispr.crispr_calculator(query_virus_dir, intermediate_dir, numThreads)
        self.crispr = src.crispr.uniGenus(self._crispr_signals, self._virus_index, self._host_index)
        # wish , if statement?
        #self.blast = src.blast.blast_calculator(query_virus_dir, self._virus_index, self._host_index, numThreads)
        
            
    def getScores(self):
        '''
        Output a prediction table
        '''
        print("Calculating prediction scores...")
        #pred = pd.DataFrame(index=virus_index, columns=host_index).fill
        if self._short:
            score = REGRESSION_COEFFICIENTS_SHORT[0]*pd.DataFrame(index=self._virus_index, columns=self._host_index).fillna(1) +  \
            REGRESSION_COEFFICIENTS_SHORT[1]*self.wish +  \
            REGRESSION_COEFFICIENTS_SHORT[2]*self.posSV +  \
            REGRESSION_COEFFICIENTS_SHORT[3]*self.negSV +  \
            REGRESSION_COEFFICIENTS_SHORT[4]*self.crispr +  \
            REGRESSION_COEFFICIENTS_SHORT[5]*self.blast
        else:
            score = REGRESSION_COEFFICIENTS[0]*pd.DataFrame(index=self._virus_index, columns=self._host_index).fillna(1) +  \
            REGRESSION_COEFFICIENTS[1]*self.s2star +  \
            REGRESSION_COEFFICIENTS[2]*self.posSV +  \
            REGRESSION_COEFFICIENTS[3]*self.negSV +  \
            REGRESSION_COEFFICIENTS[4]*self.crispr +  \
            REGRESSION_COEFFICIENTS[5]*self.blast
        self.score = 1- 1/(pd.DataFrame(score, dtype=np.float).apply(np.exp)+1)
        
    def prediction(self, topN, output_dir_pred):  
        # dict_pred = {}
        # read in taxa info:
        taxa_info = pd.read_pickle(TAXA_INFO)
        taxa_info = taxa_info.set_index('hostNCBIName')
        '''
        If the highest scores have a tie more than topN, output them all
        '''
        print("Making predictions...")
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
            pred = taxa_info.loc[topIdx]
            pred['score'] = self.score.loc[query][topIdx]
            if self._short:
                pred['WIsH_val'] = self.wish.loc[query]
            else:
                pred['s2star_val'] = self.s2star.loc[query]
            pred['posSV_val'] = self.posSV.loc[query]
            pred['negSV_val'] = self.negSV.loc[query]
            pred['crispr_val'] = self.crispr.loc[query]
            pred['blast_val'] = self.blast.loc[query]
            if self._short:
                pred['WIsH_pct'] = self.wish.loc[query].rank(pct=True, method='min').loc[topIdx]
            else:
                pred['s2star_pct'] = self.s2star.loc[query].rank(pct=True, method='min').loc[topIdx]
            if self.posSV.loc[query].max() == 0:
                pred['posSV_pct'] = ['NA'] * topNum 
            else: 
                pred['posSV_pct'] = self.posSV.loc[query].rank(pct=True, method='min').loc[topIdx]
            if self.negSV.loc[query].max() == 0:
                pred['negSV_pct'] = ['NA'] * topNum
            else:     # negative coefficient
                pred['negSV_pct'] = 1 - self.negSV.loc[query].rank(pct=True, method='max').loc[topIdx]
            if self.crispr.loc[query].max() == 0:
                pred['crispr_pct'] = ['NA'] * topNum
            else: 
                pred['crispr_pct'] = self.crispr.loc[query].rank(pct=True, method='min').loc[topIdx]
            if self.blast.loc[query].max() == 0:
                pred['blast_pct'] = ['NA'] * topNum
            else:
                pred['blast_pct'] = self.blast.loc[query].rank(pct=True, method='min').loc[topIdx]
            #dict_pred[query] = pred
            pred.to_csv(os.path.join(output_dir_pred, (query+'_prediction.csv')), float_format='%.4f')
        #return dict_pred
            
            
        

'''
Main function
'''
import argparse
import os
import sys

parser = argparse.ArgumentParser(description='VHiNet: Predict hosts given query viruses')
parser.add_argument('-q', dest='query_virus_dir',nargs=1,required=True, help='Directory containing query virus genomes with .fasta or .fa suffix')
parser.add_argument('-t',dest='num_Threads',nargs=1,type=int,default=[1], help='Number of threads (CPUs) to use in the BLAST search. Default = 1')
parser.add_argument('--short-contig',action='store_true',help='Predict hosts for short viral contigs. WIsH model files are required in this mode')
parser.add_argument('-o',dest='output_dir',nargs=1,required=True,help='Output directory')
parser.add_argument('-n',dest='topN',metavar='topN',nargs=1,type=int,default=[1], help='Number of top predictions written to the output files. All predictions will be output if there is a tie at the highest score. Default = 1')

args = parser.parse_args()
#args = parser.parse_args(['-q','test_query/','-t','8','-o','tmp','--short-contig'])

query_virus_dir = os.path.abspath(os.path.expanduser(args.query_virus_dir[0]))
if not os.path.isdir(query_virus_dir):
    sys.exit('Query directory error: no such directory')


output_dir = os.path.abspath(os.path.expanduser(args.output_dir[0]))
if not os.path.isdir(output_dir):
    sys.exit('Output directory error: no such directory')

# =============================================================================
# Prediction 
# =============================================================================
print('Loading packages...')

import pandas as pd    
from predictor import HostPredictor     

predictor = HostPredictor(query_virus_dir, args.short_contig, args.num_Threads[0])

output_dir_features = os.path.join(output_dir, 'feature_values')
try:
    os.stat(output_dir_features)
except:
    os.mkdir(output_dir_features)


if args.short_contig:
    predictor.wish.to_csv(os.path.join(output_dir_features,'feature_values_wish.csv'))
else:
    predictor.s2star.to_csv(os.path.join(output_dir_features,'feature_values_s2star.csv'))


predictor.crispr.to_csv(os.path.join(output_dir_features,'feature_values_crispr.csv'))
predictor.posSV.to_csv(os.path.join(output_dir_features,'feature_values_posSV.csv'))
predictor.negSV.to_csv(os.path.join(output_dir_features,'feature_values_negSV.csv'))


predictor.getScores()

dict_pred = predictor.prediction(args.topN[0])
for query,preds in dict_pred.items():
    preds.to_csv(os.path.join(output_dir, (query+'_prediction.csv')))

print('---- Predictions were written to ',output_dir,' ----')
    



























#!/usr/bin/env python
# =============================================================================
# Sanity check
# =============================================================================
import argparse
import os
import sys

parser = argparse.ArgumentParser(description='VirHostMatcher-Net: A network-based tool for predicting hosts given query viruses')
parser.add_argument('-q', dest='query_virus_dir',nargs=1,required=True, help='Directory containing query virus genomes with .fasta or .fa suffix')
parser.add_argument('-o',dest='output_dir',nargs=1,required=True,help='Output directory')
parser.add_argument('-t',dest='num_Threads',nargs=1,type=int,default=[1], help='Number of threads to use. Default = 1')
parser.add_argument('--short-contig',action='store_true',help='Predict hosts for short viral contigs. WIsH model files are required in this mode')
parser.add_argument('-n',dest='topN',metavar='topN',nargs=1,type=int,default=[1], help='Number of top predictions written to the output files. All predictions will be output if there is a tie in score. Default = 1')
parser.add_argument('-i',dest='intermediate_dir',nargs=1, default=['./intermediate_res'], help='Directory storing intermediate result. Default = ./intermediate_res')
parser.add_argument('-l',dest='genome_list',nargs=1, default=[None], help='Location of the file containing host NCBI genome names of interest')

args = parser.parse_args()
#args = parser.parse_args(['-q','test_query/','-t','8','-o','tmp','--short-contig'])

query_virus_dir = os.path.abspath(os.path.expanduser(args.query_virus_dir[0]))
if not os.path.isdir(query_virus_dir):
    sys.exit('Query directory error: no such directory '+ query_virus_dir)


output_dir = os.path.abspath(os.path.expanduser(args.output_dir[0]))
if not os.path.isdir(output_dir):
    sys.exit('Output directory error: no such directory '+ output_dir)

genome_list = args.genome_list[0]
if genome_list is not None:
    if not os.path.isfile(genome_list):
        sys.exit('Genome ID file error: no such file '+ genome_list)

intermediate_dir = os.path.abspath(os.path.expanduser(args.intermediate_dir[0]))
    

# =============================================================================
# Prediction 
# =============================================================================
print('Loading packages...')

import pandas as pd    
from predictor import HostPredictor     

predictor = HostPredictor(query_virus_dir, args.short_contig, intermediate_dir, genome_list, args.num_Threads[0])

output_dir_features = os.path.join(output_dir, 'feature_values')
try:
    os.stat(output_dir_features)
except:
    os.mkdir(output_dir_features)

output_dir_pred = os.path.join(output_dir, 'predictions')
try:
    os.stat(output_dir_pred)
except:
    os.mkdir(output_dir_pred)

    
if args.short_contig:
    predictor.wish.to_csv(os.path.join(output_dir_features,'feature_values_wish.csv'), float_format='%.4f')
else:
    predictor.s2star.to_csv(os.path.join(output_dir_features,'feature_values_s2star.csv'), float_format='%.4f')

print('Writing feature scores to {}...'.format(output_dir_features))
predictor.crispr.to_csv(os.path.join(output_dir_features,'feature_values_crispr.csv'), float_format='%.4f')
predictor.posSV.to_csv(os.path.join(output_dir_features,'feature_values_posSV.csv'), float_format='%.4f')
predictor.negSV.to_csv(os.path.join(output_dir_features,'feature_values_negSV.csv'), float_format='%.4f')
# predictor.blast.to_csv(os.path.join(output_dir_features,'feature_values_blast.csv'), float_format='%.5f')

predictor.getScores()
predictor.prediction(args.topN[0], output_dir_pred)

## write predictions
# for query,preds in dict_pred.items():
#     preds.to_csv(os.path.join(output_dir_pred, (query+'_prediction.csv')))

print('---- Predictions are written to ',output_dir,' ----')
    



























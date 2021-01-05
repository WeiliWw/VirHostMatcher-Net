#!/usr/bin/env python
# =============================================================================
# Sanity check
# =============================================================================
import os
import sys
import src.args

args = src.args.parse_arguments()
#args = parser.parse_args(['-q','test_query/','-t','8','-o','tmp','--short-contig'])

query_virus_dir = os.path.abspath(os.path.expanduser(args.query_virus_dir[0]))
if not os.path.isdir(query_virus_dir):
    sys.exit('Query directory error: no such directory '+ query_virus_dir)


output_dir = os.path.abspath(os.path.expanduser(args.output_dir[0]))
if not os.path.isdir(output_dir):
    print("Specified output directory does not exist!\n"
          "Creating {}".format(output_dir)
          )
    os.makedirs(output_dir)

genome_list = args.genome_list[0]
if genome_list is not None:
    if not os.path.isfile(genome_list):
        sys.exit('Genome ID file error: no such file '+ genome_list)

intermediate_dir = os.path.abspath(os.path.expanduser(args.intermediate_dir[0]))


# Basic checks on the data directory
data_dir = args.data_dir
if not os.path.isdir(data_dir):
    sys.exit(
        "Please provide a valid path to a directory with all required data"
    )

for subname in ['tables', 'crispr_db_prefix', 'host_wish_model']:
    subdir = os.path.join(data_dir, subname)
    if not os.path.isdir(subdir):
        sys.exit(
            "\nError: The data directory you provided is malformatted.\n"
            "More info: "
            "https://github.com/WeiliWw/VirHostMatcher-Net#downloading\n"
        )

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

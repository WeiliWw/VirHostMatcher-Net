import os
import sys

from .Variables import set_data_dir
from .config import RuntimeConfig


def _ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def validate_runtime_config(config: RuntimeConfig):
    if not os.path.isdir(config.query_virus_dir):
        sys.exit('Query directory error: no such directory ' + config.query_virus_dir)

    if config.genome_list is not None and not os.path.isfile(config.genome_list):
        sys.exit('Genome ID file error: no such file ' + config.genome_list)

    if not os.path.isdir(config.data_dir):
        sys.exit("Please provide a valid path to a directory with all required data")

    for subname in ['tables', 'crispr_db_prefix', 'host_wish_model']:
        subdir = os.path.join(config.data_dir, subname)
        if not os.path.isdir(subdir):
            sys.exit(
                "\nError: The data directory you provided is malformatted.\n"
                "More info: "
                "https://github.com/WeiliWw/VirHostMatcher-Net#download\n"
            )


def run_prediction(config: RuntimeConfig):
    validate_runtime_config(config)

    if not os.path.isdir(config.output_dir):
        print(
            "Specified output directory does not exist!\n"
            "Creating {}".format(config.output_dir)
        )
        _ensure_dir(config.output_dir)

    set_data_dir(config.data_dir)

    print('Loading packages...')
    from predictor import HostPredictor

    predictor = HostPredictor(
        config.query_virus_dir,
        config.short_contig,
        config.intermediate_dir,
        config.genome_list,
        config.num_threads,
    )

    output_dir_features = os.path.join(config.output_dir, 'feature_values')
    output_dir_pred = os.path.join(config.output_dir, 'predictions')
    _ensure_dir(output_dir_features)
    _ensure_dir(output_dir_pred)

    if config.short_contig:
        predictor.wish.to_csv(
            os.path.join(output_dir_features, 'feature_values_wish.csv'),
            float_format='%.4f',
        )
    else:
        predictor.s2star.to_csv(
            os.path.join(output_dir_features, 'feature_values_s2star.csv'),
            float_format='%.4f',
        )

    print('Writing feature scores to {}...'.format(output_dir_features))
    predictor.crispr.to_csv(
        os.path.join(output_dir_features, 'feature_values_crispr.csv'),
        float_format='%.4f',
    )
    predictor.posSV.to_csv(
        os.path.join(output_dir_features, 'feature_values_posSV.csv'),
        float_format='%.4f',
    )
    predictor.negSV.to_csv(
        os.path.join(output_dir_features, 'feature_values_negSV.csv'),
        float_format='%.4f',
    )

    predictor.getScores()
    predictor.prediction(config.topN, output_dir_pred)
    print('---- Predictions are written to ', config.output_dir, ' ----')

# Make arguments accessible from modules
# see https://stackoverflow.com/a/41324084
import argparse
import os

def parse_arguments():
    # Define default data dir to be the same as before
    up_one = os.path.dirname(os.path.dirname((os.path.abspath(__file__))))
    # Append data to it
    DEFAULT_DATA_DIR = os.path.join(up_one, 'data')

    parser = argparse.ArgumentParser(
        description="VirHostMatcher-Net: A network-based tool for predicting "
        "hosts given query viruses",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    optionalArgs = parser._action_groups.pop()
    optionalArgs.title= "Optional Arguments"

    requiredArgs = parser.add_argument_group("Required Arguments")
    requiredArgs.add_argument('-q',
                        dest='query_virus_dir',
                        required=True,
                        type=str,
                        help="Directory containing query virus genomes with "
                        ".fasta, .fa or .fna  suffixes. One genome per file."
                        )
    requiredArgs.add_argument('-o',
                        dest='output_dir',
                        type=str,
                        required=True,
                        help="Output directory. It is created if it doesn't "
                        "exist"
                        )
    optionalArgs.add_argument('-t',
                        dest='num_Threads',
                        required=False,
                        type=int,
                        default=1,
                        help='Number of threads to use.'
                        )
    optionalArgs.add_argument('--short-contig',
                        action='store_true',
                        required=False,
                        help="Predict hosts for short viral contigs. "
                        "(WIsH model files are required in this mode"
                        )
    optionalArgs.add_argument('-n',
                        dest='topN',
                        metavar='topN',
                        required=False,
                        type=int,
                        default=1,
                        help="Number of top predictions written to the output "
                        "files. All predictions will be output if there is a tie "
                        "in score"
                        )
    optionalArgs.add_argument('-i',
                        dest='intermediate_dir',
                        default='./intermediate_res',
                        type=str,
                        required=False,
                        help="Directory storing intermediate result."
                        )
    optionalArgs.add_argument('-l',
                        dest='genome_list',
                        required=False,
                        default=None,
                        help="Location of the file containing host "
                        "NCBI genome names of interest"
                        )
    optionalArgs.add_argument('-d',
                        dest='data_dir',
                        type=str,
                        required=False,
                        default=DEFAULT_DATA_DIR,
                        help="Directory where models, blast and CRISPRs are "
                        "located."
                        )

    parser._action_groups.append(optionalArgs)

    return parser.parse_args()

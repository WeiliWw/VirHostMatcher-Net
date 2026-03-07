'''
Preset variables
'''
import os

# Default data directory (no import-time CLI parsing).
_DEFAULT_DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'data',
)
_path = _DEFAULT_DATA_DIR

# intercept, s2star/wish, posSV, negSV, crispr
REGRESSION_COEFFICIENTS = [-1.85696213,  17.13324415, 4.27680508, -27.73515517, 0.12533108]
REGRESSION_COEFFICIENTS_SHORT = [39.58978537,  27.69902541, 6.40628182, -15.18968327, 0.20333482]
# INTERMEDIATE_RESULT = os.path.join(_path, 'intermediate_res/')
DB_HOST_PREFIX = os.path.join(_path, 'host_db_prefix/hostGenomes.fna')
DB_HOST_CRISPR_PREFIX = os.path.join(_path, 'crispr_db_prefix/allCRISPRs')
HASH_TABLE = os.path.join(_path, 'tables/genome2header.pkl')
WISH_HOST_MODELS =  os.path.join(_path, 'host_wish_model/')
TAXA_INFO = os.path.join(_path, 'tables/hostTaxa.pkl')
PRED_THRE = os.path.join(_path, 'tables/predThre.csv')
# TABLES: intermediate s2 values and the interaction 0-1 matrix
# TABLES = os.path.join(_path, 'data/tables/tables.h5')

TABLE_HOST = os.path.join(_path, 'tables/host62493_s2_mat.pkl')
TABLE_BENCH = os.path.join(_path, 'tables/bench_s2_mat.pkl')
TABLE_INTER = os.path.join(_path, 'tables/df_interaction.csv')

# D2STAR_HASH = os.path.join(_path, 'data/hash/')
# S2STAR_BENCH_HOST = os.path.join(_path, 'data/tables/relationMat_352by31986.csv')
# PSEUDO_HOST = os.path.join(_path, 'data/pseudo_host/')
# PSEUDO_VIRUS = os.path.join(_path, 'data/pseudo_virus/')

# try:
#     os.stat(INTERMEDIATE_RESULT)
# except:
#     os.mkdir(INTERMEDIATE_RESULT)


def set_data_dir(data_dir):
    """
    Configure runtime data directory explicitly.
    Must be called before importing modules that consume these constants.
    """
    global _path
    global DB_HOST_PREFIX, DB_HOST_CRISPR_PREFIX, HASH_TABLE, WISH_HOST_MODELS
    global TAXA_INFO, PRED_THRE, TABLE_HOST, TABLE_BENCH, TABLE_INTER

    _path = os.path.abspath(os.path.expanduser(data_dir))
    DB_HOST_PREFIX = os.path.join(_path, 'host_db_prefix/hostGenomes.fna')
    DB_HOST_CRISPR_PREFIX = os.path.join(_path, 'crispr_db_prefix/allCRISPRs')
    HASH_TABLE = os.path.join(_path, 'tables/genome2header.pkl')
    WISH_HOST_MODELS = os.path.join(_path, 'host_wish_model/')
    TAXA_INFO = os.path.join(_path, 'tables/hostTaxa.pkl')
    PRED_THRE = os.path.join(_path, 'tables/predThre.csv')
    TABLE_HOST = os.path.join(_path, 'tables/host62493_s2_mat.pkl')
    TABLE_BENCH = os.path.join(_path, 'tables/bench_s2_mat.pkl')
    TABLE_INTER = os.path.join(_path, 'tables/df_interaction.csv')

'''
Preset variables
'''
import os
_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# intercept, s2star/wish, posSV, negSV, crispr, blast
# REGRESSION_COEFFICIENTS = [-0.648915, 14.756166, 3.401181,
#                            -29.188570, 0.324167, 6334.783871]
REGRESSION_COEFFICIENTS = [-1.85696213,  17.13324415, -27.73515517, 4.27680508, 0.12533108]
# REGRESSION_COEFFICIENTS_SHORT = [28.311335, 19.282953, 5.485007,
#                             -17.328643, 0.412549, 4462.823622]
REGRESSION_COEFFICIENTS_SHORT = [39.58978537,  27.69902541, -15.18968327, 6.40628182, 0.20333482]
# INTERMEDIATE_RESULT = os.path.join(_path, 'intermediate_res/')
# DB_HOST_PREFIX = os.path.join(_path, 'data/host_db_prefix/hostGenomes.fna')
DB_HOST_CRISPR_PREFIX = os.path.join(_path, 'data/crispr_db_prefix/allCRISPRs')
HASH_TABLE = os.path.join(_path, 'data/tables/genome2header.pkl')
WISH_HOST_MODELS =  os.path.join(_path, 'data/host_wish_model/')
TAXA_INFO = os.path.join(_path, 'data/tables/hostTaxa.pkl')
PRED_THRE = os.path.join(_path, 'data/tables/predThre.csv')
# TABLES: intermediate s2 values and the interaction 0-1 matrix
# TABLES = os.path.join(_path, 'data/tables/tables.h5')

TABLE_HOST = os.path.join(_path, 'data/tables/host62493_s2_mat.pkl')
TABLE_BENCH = os.path.join(_path, 'data/tables/bench_s2_mat.pkl')
TABLE_INTER = os.path.join(_path, 'data/tables/df_interaction.csv')

# D2STAR_HASH = os.path.join(_path, 'data/hash/')
# S2STAR_BENCH_HOST = os.path.join(_path, 'data/tables/relationMat_352by31986.csv')
# PSEUDO_HOST = os.path.join(_path, 'data/pseudo_host/')
# PSEUDO_VIRUS = os.path.join(_path, 'data/pseudo_virus/')

# try:
#     os.stat(INTERMEDIATE_RESULT)
# except:
#     os.mkdir(INTERMEDIATE_RESULT)

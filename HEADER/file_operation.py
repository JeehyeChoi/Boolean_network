# file_operations.py

import os
import pandas as pd
import numpy as np

def import_truth_table(fnorigin, KEY):
    fn_tt = "DATA_{:s}/{:s}_TT/{:s}.csv".format("truth", fnorigin, KEY)
    if os.path.exists(fn_tt):
        df_tt = pd.read_csv(fn_tt, header=0, sep=", ", dtype=np.float64, engine='python')
    else:
        df_tt = pd.DataFrame(columns=['aa'])
    return df_tt

def import_key(fnorigin):
    fn_skey = "DATA_{:s}/{:s}_TT/SPECIES_KEY.csv".format("truth", fnorigin)
    if os.path.exists(fn_skey):
        df_skey = pd.read_csv(fn_skey, header=None, sep="\t")
        df_skey.columns = ['Value_Column', 'Key_Column']
        df_skey = df_skey.set_index('Key_Column')
    else:
        df_skey = pd.DataFrame()
    return df_skey

def import_adj(fnorigin):
    fn = "DATA_{:s}/{:s}_MATRIX_pm.csv".format("matrix", fnorigin)
    df_matrix = pd.read_csv(fn, header=0, sep=",")
    adj_header = df_matrix.columns
    df_matrix.index = adj_header
    adj_matrix = df_matrix
    return adj_header, adj_matrix


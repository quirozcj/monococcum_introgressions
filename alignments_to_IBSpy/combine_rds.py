import argparse
import pandas as pd
import numpy as np
from itertools import count
import pyreadr
from warnings import filterwarnings
filterwarnings('ignore')
from datetime import datetime

### Script to rename and compile all chromosome and references divided by query sample for coverage analysis ###

startTime = datetime.now()

def map_substring(s, dict_map):
    for key in dict_map.keys():
        if key in s:
            return dict_map[key]
    return s

def rename_header(df, names_dic):
    df_r = df.rename(columns={'rid': 'chr', 'r_length':'aln_length','perc_id':'aln_perc_id','rs': 'r_start', 're': 'r_end', 'qs': 'q_start', 'qe': 'q_end', 'qid': 'query'}).reset_index(drop=True)
    df_ord = df_r[['chr', 'r_start', 'r_end', 'aln_length', 'query', 'q_start', 'q_end', 'aln_perc_id', 'perc_id_factor', 'r_mid', 'q_mid','error', 'strand','comparison']]
    
    ref_col = df_ord['chr'].apply(lambda x: map_substring(x, names_dic))
    df_ord.insert(loc=3, column='reference', value=ref_col)
    
    query_col = df_ord['query'].apply(lambda x: map_substring(x, names_dic))
    df_ord.insert(loc=5, column='query_aln', value=query_col)
    df_ord[['r_start', 'r_end']] = df_ord[['r_start', 'r_end']].astype(int)
    return df_ord

def exctract_query(df, query):
    df_q = df[df['query_aln'] == query]
    return df_q

def combine_db(rds_files):
    rds_df = pd.read_csv(rds_files, delimiter='\t')
    dfs = []
    for index, row in rds_df.iterrows():
        in_db = pyreadr.read_r(row['path']+row['file'])[None]
        dfs.append(in_db)
    dfs_out = pd.concat(dfs, axis=0)
    return dfs_out

def create_dic(pfx_file):
    with open(pfx_file) as f:
     rows = (line.rstrip().split('\t') for line in f)
     dic = {row[0]:row[1] for row in rows}
     return dic

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--rds_files')
    parser.add_argument('-q', '--query_aln')
    parser.add_argument('-p', '--pfx_file')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    names_dic = create_dic(args.pfx_file)
    dfs_combined = combine_db(args.rds_files)
    df_renamed = rename_header(dfs_combined, names_dic)
    df_query = exctract_query(df_renamed, args.query_aln)
    df_query.to_csv(args.output, sep='\t', index=False)
    print(datetime.now() - startTime)

if __name__ == '__main__':
    main()
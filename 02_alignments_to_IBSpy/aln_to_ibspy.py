import argparse
import pandas as pd
import numpy as np
from itertools import count
from datetime import datetime
startTime = datetime.now()

## functions start ###
def map_substring(s, dict_map):
    for key in dict_map.keys():
        if key in s:
            return dict_map[key]
    return s

def get_chr_len(file, reference, chromosome):
    len_df = pd.read_csv(file, delimiter='\t')
    ref_df = len_df[(len_df['reference'] == reference.lower()) & len_df['chr'].str.contains(chromosome)].iloc[0][3]
    return ref_df

def exctract_rf(df, reference):
    df_r = df[df['reference'] == reference.lower()]
    return df_r

def bin_data(df, chr_id, chr_len, w_size):
    df = df[df['chr'].str.contains(chr_id)].copy()
    w_pos=1
    dfs = []
    for i in count(w_pos, w_size):
        if i > chr_len:
            break
        windows_df = df[(df['r_mid'] >= w_pos) & (df['r_mid'] <= w_pos + w_size)].copy()
        if not windows_df.empty:
            windows_df.insert(1, 'start', w_pos)
        else:
            windows_df = pd.DataFrame(columns=df.columns, index=range(1))

            windows_df['chr'] = df['chr'].iloc[0]
            windows_df.insert(1, 'start', w_pos)
            windows_df.fillna(0, inplace=True)
        dfs.append(windows_df)
        w_pos += w_size # add sliding windows
    dfs_cat = pd.concat(dfs, axis=0)
    dfs_cat.insert(2, 'end', dfs_cat['start'] + w_size - 1)
    dfs_cat['query_aln'] = df['query_aln'].iloc[0]
    dfs_cat['reference'] = df['reference'].iloc[0]
    dfs_cat['comparison'] = df['comparison'].iloc[0]
    return dfs_cat

def stats_alignment(df):
    df_g = df.groupby(['chr','start','end','query_aln'])
    df_g_mid = pd.DataFrame(df_g[['aln_perc_id','aln_length']].median()).reset_index()
    return df_g_mid

def count_by_windows(df, chr_id, chr_len, w_size):
    df_chr = df[df['seqname'].str.contains(chr_id)].copy()
    w_pos=0    
    dfs = []
    for i in count(w_pos, w_size):
        if i > chr_len:
            break
        by_windows_df = df_chr[(df_chr['end'] > w_pos) & (df_chr['end'] <= w_pos + w_size)].copy()
        by_windows_df['start'] = w_pos + 1
        dfs.append(by_windows_df)
        w_pos += w_size
    dfs_cat = pd.concat(dfs, axis=0)
    dfs_cat['end'] = dfs_cat['start'] + w_size - 1
    return dfs_cat

def stats_ibs(df):
    df_g = df.groupby(['seqname','start','end']).sum().reset_index()
    df_g['kmr_perc_id'] = df_g['observed_kmers']/df_g['total_kmers']*100
    df_g['50kb_variations'] = df_g['variations']/(df_g.iloc[0][2]/50000)
    return df_g

def rename_bed(df, dict_map):
    ref_col = df[0].apply(lambda x: map_substring(x, dict_map))
    df.insert(loc=3, column='reference', value=ref_col)
    return df

def filter_bed(df, reference, chromosome):
    df_r = df[(df['reference'] == reference.lower()) & (df[0].str.contains(chromosome))].reset_index(drop=True)
    return df_r


def combine_dfs(df1, df2, query_ibs):
    df1['variations'] = df2['variations']
    df1['kmr_perc_id'] = df2['kmr_perc_id']
    df1['50kb_variations'] = df2['50kb_variations']
    df1.insert(9, 'query_ibs', query_ibs)
    return df1

def create_dic(pfx_file):
    with open(pfx_file) as f:
     rows = (line.rstrip().split('\t') for line in f)
     dic = {row[0]:row[1] for row in rows}
     return dic
### end of functions ###


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--ibs_path')
    parser.add_argument('-a', '--aln_file')
    parser.add_argument('-v', '--cov_file')
    parser.add_argument('-r', '--reference')
    parser.add_argument('-b', '--query_ibs')
    parser.add_argument('-w', '--window', type=int)
    parser.add_argument('-c', '--chromosomes', nargs='+')
    parser.add_argument('-l', '--chr_lenghts', type=str)
    parser.add_argument('-p', '--pfx_file')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    reference = args.reference
    w_size = args.window
    query_ibs = args.query_ibs

    names_dic = create_dic(args.pfx_file)
    df_renamed = pd.read_csv(args.aln_file, delimiter='\t')
    coverage_f = pd.read_csv(args.cov_file, delimiter='\t', header=None)
    coverage_rn = rename_bed(coverage_f, names_dic)


    dfs = []
    for chromosome in args.chromosomes:
        chr_len = get_chr_len(args.chr_lenghts, reference, chromosome)
        if reference == 'chinesespring':
            aln_reference = 'chinese'
            rq_df = exctract_rf(df_renamed, aln_reference)
            coverage_df = filter_bed(coverage_rn, aln_reference, chromosome)
        else:
            rq_df = exctract_rf(df_renamed, reference)
            coverage_df = filter_bed(coverage_rn, reference, chromosome)

        dfs_cat = bin_data(rq_df, chromosome, chr_len, w_size)
        aln_stats = stats_alignment(dfs_cat)
        aln_stats[['coverage','coverage_prc','depth']] = coverage_df[[4,6,3]]

        ibs_file = pd.read_csv(args.ibs_path+f'/{query_ibs}_{reference}_50000.tsv.gz', delimiter='\t')
        ibspy_df = count_by_windows(ibs_file, chromosome, chr_len, w_size)
        ibs_stats = stats_ibs(ibspy_df)
        out_df = combine_dfs(aln_stats, ibs_stats, query_ibs)
        dfs.append(out_df)

    dfs_combined = pd.concat(dfs, axis=0)
    print(dfs_combined)
    dfs_combined.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()
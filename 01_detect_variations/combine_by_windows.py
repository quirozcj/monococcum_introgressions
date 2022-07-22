import argparse
import pandas as pd
from datetime import datetime

### Scritp to combine multiple IBSpy output files by window & reference  ###

startTime = datetime.now()

def combine_db(samples_id, reference, score):
    samples_id = pd.read_csv(samples_id, delimiter='\t')
    samples_by_reference = samples_id[samples_id['reference'] == reference]

    samples_combined_db = pd.DataFrame()
    for index, row in samples_by_reference.iterrows():
        in_db = pd.read_csv(row['path']+row['file'], delimiter='\t')
        sample_name = row['query']
        if score == 'observed_kmers':
            samples_combined_db[sample_name] = in_db['observed_kmers']/in_db['total_kmers']
        else:
            samples_combined_db[sample_name] = in_db['variations']
            
    row_names = in_db[['seqname', 'start', 'end']]
    samples_combined_db = pd.concat([row_names, samples_combined_db], axis=1)
    return samples_combined_db

def count_by_windows(combined_samples, window_size, score):
    in_db = combined_samples
    window_size = window_size
    chrLen = in_db['end'].max()
    
    w_pos = 0
    db_byChr = pd.DataFrame()
    while w_pos <= chrLen:
        by_windows_df = in_db[(in_db['end'] > w_pos) & (in_db['end'] <= w_pos + window_size)]
        by_windows_df = by_windows_df.drop(['start','end'], axis=1)
        by_windows_df['start'] = w_pos + 1
        by_windows_df['end'] = w_pos + window_size
        w_pos += window_size
        db_byChr = db_byChr.append(by_windows_df)  
    if score == 'observed_kmers':
        by_windows_db = db_byChr.groupby(['seqname','start','end']).mean().reset_index()
    else:
        by_windows_db = db_byChr.groupby(['seqname','start','end']).sum().reset_index()
    return by_windows_db

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--samples_id')
    parser.add_argument('-w', '--window_size', type=int)
    parser.add_argument('-r', '--reference')
    parser.add_argument('-s', '--score')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    db_combined = combine_db(args.samples_id, args.reference, args.score)
    by_windowd_db = count_by_windows(db_combined, args.window_size, args.score)
    by_windowd_db.to_csv(args.output, sep='\t', index=False)
    print(datetime.now() - startTime)
if __name__ == '__main__':
    main()
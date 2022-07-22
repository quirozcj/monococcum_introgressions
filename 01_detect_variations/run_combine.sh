#!/bin/bash
#SBATCH --partition=jic-long
#SBATCH --nodes=1
#SBATCH --cpus=1
#SBATCH --mem 20G
#SBATCH -o /log/combine.%N.%j.out # STDOUT
#SBATCH -e /log/combine.%N.%j.err # STDERR
#SBATCH --job-name=combine

reference='reference_name'
window=50000
score='variations'

out_dir=../ibspy_combined
mkdir -p $out_dir

python3 combine_by_windows.py \
-i metadata.tsv \
-r ${reference} \
-w ${window} \
-s ${score} \
-o $out_dir/${reference}_combined_queries_${window}w.tsv.gz
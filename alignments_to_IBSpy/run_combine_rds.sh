#!/bin/bash
#SBATCH --partition=jic-short
#SBATCH --nodes=1
#SBATCH --cpus=1
#SBATCH --mem 10G
#SBATCH -o /log/combine_rds.%N.%j.out # STDOUT
#SBATCH -e /log/combine_rds.%N.%j.err # STDERR
#SBATCH --job-name=combine_rds
#SBATCH --array=0-14

### Run:Script to rename and compile all chromosome and references divided by query sample for coverage analysis ###

i=$SLURM_ARRAY_TASK_ID

declare -a queries_aln=(\
	'arinalrfor' \
	'symattis' \
	'julius' \
	'mace' \
	'norin61' \
	'stanley' \
	'landmark' \
	'lancer' \
	'jagger' \
	'chinese' \
	'TA299' \
	'TA10622' \
	'urartu' \
	'svevo' \
	'zavitan' \
	)

out_dir=../09_watseq/15_alignments_to_ibspy/mummer_rds_combined/
mkdir -p ${out_dir}

singularity exec ~/tmp/quirozc/python3.img python3 combine_rds.py \
-f rds_files.tsv \
-q ${query_aln} \
-p query_prefix.tsv \
-o ${out_dir}/whole_genome_aln_all_references_vs_${query_aln}.tsv
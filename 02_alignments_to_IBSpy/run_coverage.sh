#!/bin/bash
#SBATCH --partition=jic-short
#SBATCH --nodes=1
#SBATCH --cpus=1
#SBATCH --mem 10G
#SBATCH -o /log/coverage.%N.%j.out # STDOUT
#SBATCH -e /log/coverage.%N.%j.err # STDERR
#SBATCH --job-name=coverage
#SBATCH --array=0-14

### Scritp to obtain the coverage by query sample of against all pangenome assemblies by window size ###

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

window=500000
rds_dir=../mummer_rds_combined

source bedtools-2.28.0

bedtools makewindows \
-g chr_lenghts_bed.tsv \
-w ${window} | \
bedtools coverage \
-a - \
-b ${rds_dir}/whole_genome_aln_all_references_vs_${query_aln}.tsv \
> ${rds_dir}/whole_genome_aln_all_references_vs_${query_aln}_${window}w_coverage.tsv
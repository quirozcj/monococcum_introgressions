#!/bin/bash
#SBATCH --partition=jic-short
#SBATCH --nodes=1
#SBATCH --cpus=1
#SBATCH --mem 10G
#SBATCH -o /log/aln_to_ibs.%N.%j.out # STDOUT
#SBATCH -e /log/aln_to_ibs.%N.%j.err # STDERR
#SBATCH --job-name=aln_to_ibs
#SBATCH --array=0-9

### Scritp to merge alignments coverage tables with IBSpy variaitons data by window ###

i=$SLURM_ARRAY_TASK_ID

declare -a references=(\
	'arinaLrFor' \
	'symattis' \
	'julius' \
	'mace' \
	'norin61' \
	'stanley' \
	'landmark' \
	'lancer' \
	'chinesespring' \
	'jagger' \
	)
reference_i=$(($i%10))
reference=${references[$reference_i]}

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
	'zavitan' \
	'TA299' \
	'TA10622' \
	'urartu' \
	'svevo' \
	)

declare -a queries_ibs=(\
	'arina-pg' \
	'mattis-pg' \
	'julius-pg' \
	'mace-pg' \
	'norin61-pg' \
	'stanley-pg' \
	'landmark-pg' \
	'lancer-pg' \
	'jagger-pg' \
	'chinese-pg' \
	'zavitan-pg' \
	'TA299-pg' \
	'TA10622-pg' \
	'urartu-pg' \
	'svevo-pg' \
	)
query_i=$(($i/10))
query_aln=${queries_aln[$query_i]}
query_ibs=${queries_ibs[$query_i]}

window=500000
declare -a chromosomes=(
	'chr1A' \
	'chr2A' \
	'chr3A' \
	'chr4A' \
	'chr5A' \
	'chr6A' \
	'chr7A' \
	)


ibs_path=../IBSpy_output/${reference}
aln_file=../mummer_rds_combined
out_dir=../aln_ibs_combined/${window}
mkdir -p ${out_dir}

if [ "$reference" == "$query_aln" ]
then
	echo "Reference and Query are the same"
else

	singularity exec ~/tmp/quirozc/python3.img python3 aln_to_ibspy.py \
	-i ${ibs_path} \
	-a ${aln_file}/whole_genome_aln_all_references_vs_${query_aln}.tsv \
	-v ${aln_file}/whole_genome_aln_all_references_vs_${query_aln}_${window}w_coverage.tsv \
	-r ${reference} \
	-b ${query_ibs} \
	-w ${window} \
	-c ${chromosomes[@]} \
	-p query_prefix.tsv \
	-l chr_lengths.tsv \
	-o $out_dir/${reference}_vs_${query_aln}_vs_${query_ibs}_${window}_aln_to_ibspy.tsv
fi
#!/bin/bash
#SBATCH --partition=jic-long
#SBATCH --nodes=1
#SBATCH --cpus=1
#SBATCH --mem 40G
#SBATCH -o /log/IBSpy.%N.%j.out
#SBATCH -e /log/IBSpy.%N.%j.err
#SBATCH --array=0-32
#SBATCH -J IBSpy

i=$SLURM_ARRAY_TASK_ID

source IBSpy-0.3.1

function log_line() {
	echo $(date) "$1" >&2
}

declare -a references=(\
	"arinaLrFor" \
	"chinesespring" \
	"jagger" \
	"julius" \
	"lancer" \
	"landmark" \
	"mace" \
	"norin61" \
	"stanley" \
	"symattis" \
	"spelta" \
	)

declare -a databases=(\
	"/path/accesion1" \
	"/path/accesion2" \
	"/path/accesion3" \
	)

declare -a db_names=(\
	"accesion1" \
	"accesion2" \
	"accesion3"
	)

ref_folder="assemblies"
k_folder="path/kmers_folder/"
cd $k_folder

reference=$(($i%11))
reference=${references[$reference]}

db_i=$(($i/11))
db=${databases[$db_i]}
name=${db_names[$db_i]}

ref=${ref_folder}/${reference}.fa
window_size=50000
log_line "Running from: $PWD"
log_line "Reference: $reference"	
log_line "db: $db"
log_line "name: $name"


out_dir=${k_folder}/${reference}
mkdir -p $out_dir

IBSpy \
--kmer_size 31 \
--window_size ${window_size} \
--reference ${ref} \
--database $db \
--database_format kmc3 \
--output ${out_dir}/${name}_${reference}_${window_size}.tsv \
--compress

log_line "DONE IBSpy"
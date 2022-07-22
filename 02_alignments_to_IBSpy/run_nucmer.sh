#!/bin/bash
#SBATCH --partition=nbi-long,jic-long,RG-Cristobal-Uauy
#SBATCH --nodes=1
#SBATCH --cpus=1
#SBATCH --mem 30G
#SBATCH -o /log/mummer.%N.%j.out # STDOUT
#SBATCH -e /log/mummer.%N.%j.err # STDERR
#SBATCH --job-name=mummer
#SBATCH --array=0-54


### Run whole geneome assembly mummer alignment by chromosome & filter by 20 Kbp length ###
source mummer-3.23

chromosome="1A"
i=$SLURM_ARRAY_TASK_ID
x=$(($i % 11))
q_x=$(($i / 11))

declare -a names=(\
	"arinalrfor" \
	"chinese" \
	"jagger" \
	"julius" \
	"lancer" \
	"landmark" \
	"mace" \
	"norin61" \
	"stanley" \
	"sy_mattis" \
	"spelta" \
	)

declare -a queries=(\
	"zavitan" \
	"urartu" \
	"TA299" \
	"TA10622" \
	"tibetan" \
	)

reference=${names[$x]}
query=${queries[$q_x]}

ref_dir=../fasta/${chromosome}
query_dir=../fasta/chr${chromosome}
our_dir=../aln/$reference/chr${chromosome}
mkdir -p $our_dir


if [ "$reference" == "$query" ]
then
	echo "Reference and Query are the same"
else
	#run nucmer
	nucmer --mum --delta \
	$ref_dir/${reference}.chr${chromosome}.fa \
	$query_dir/${query}.chr${chromosome}.fa \
	--prefix $our_dir/${reference}_v_${query}.chr${chromosome} \

	#filter for alignments of at least 20kb and rq
	delta-filter -l 20000 \
	-r \
	-q $our_dir/${reference}_v_${query}.chr${chromosome}.delta \
	> $our_dir/${reference}_v_${query}.chr${chromosome}_filtered_L20Kb_rq.delta
fi

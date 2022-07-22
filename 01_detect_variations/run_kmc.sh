#!/bin/bash
#SBATCH --job-name=kmc
#SBATCH --partition=jic-short
#SBATCH --nodes=1
#SBATCH --cpus=5
#SBATCH --mem 30G

source kmc-3.0.1

#### genome assemblies ###
assembly_id=julius
in_dir=../${accesion_id}
out_dir=../out/${accesion_id}
mkdir -p $out_dir

kmc \
-k31 \
-ci1 \
-m30 \
-t5 \
-fm \
$in_dir/${assembly_id}.fa ${out_dir}/${assembly_id} ${out_dir}

#### For raw reads ###
accesion_id='accesion1'
in_dir=../${accesion_id}
out_dir=../out/${accesion_id}
mkdir -p $out_dir

kmc \
-k31 \
-ci1 \
-m30 \
-t5 \
-fq \
<(ls -d $in_dir/*.gz) ${out_dir}/${accesion_id} ${out_dir}
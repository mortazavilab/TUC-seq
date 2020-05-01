#!/bin/bash
#$ -ckpt restart
#$ -q sam
#$ -pe openmp 16

sample_dir=$1
sample_name=$2

repo_path="/dfs3/samlab/sorenar/TUC-seq"
data_path=${repo_path}'data/'${sample_dir}
result_path=${repo_path}'results/'${sample_dir}
transcriptome_ref=${repo_path}'ref/GRCh38.p12/transcriptome'


STAR --runThreadN 16 --genomeDir ${transcriptome_ref} --readFilesIn ${result_path}${sample_name}_R1.fastq --outFileNamePrefix ${result_path}${sample_name} --outFilterMismatchNmax 15 --outFilterMismatchNoverReadLmax 0.07 --outFilterMultimapNmax 10 --outSAMunmapped None --outSAMattributes MD NM --alignIntronMax 10 --alignIntronMin 20 --alignMatesGapMax 1000000 --outSAMtype BAM SortedByCoordinate



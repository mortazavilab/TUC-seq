#!/bin/bash
#$ -ckpt restart
#$ -q sam
#$ -pe openmp 1


sample_dir=$1
sample_name=$2

repo_path="/dfs3/samlab/sorenar/TUC-seq"
result_path=${repo_path}'results/'${sample_dir}
annotation_file=${repo_path}'/ref/GRCh38_transcriptome.fa'



samtools view -H ${result_path}${sample_name}_Aligned.sortedByCoord.out.bam > ${result_path}${sample_name}_header.txt
samtools view -q 2 -f 0 ${result_path}${sample_name}_Aligned.sortedByCoord.out.bam >  ${result_path}${sample_name}_filtered_1.sam
samtools view -q 2 -f 16 ${result_path}${sample_name}_Aligned.sortedByCoord.out.bam >  ${result_path}${sample_name}_filtered_2.sam
cat ${result_path}${sample_name}_header.txt ${result_path}${sample_name}_filtered_*.sam > ${result_path}${sample_name}_filtered.sam

python substitution_annotator.py -i ${result_path}${sample_name}_filtered.sam -o ${result_path}${sample_name}

for TC in 0 1 2 3
do

  samtools sort -n ${result_path}${sample_name}_TC_${TC}.sam > ${result_path}${sample_name}_${TC}_TC_sorted.sam
  express --no-bias-correct -o ${result_path} ${annotation_file} ${result_path}${sample_name}_${TC}_TC_sorted.sam
  mv ${result_path}/results.xprs ${result_path}${sample_name}_${TC}_TC_counts.txt

  samtools view -bh ${result_path}${sample_name}_TC_${TC}.sam > ${result_path}${sample_name}_${TC}_TC.bam
  samtools sort ${result_path}${sample_name}_${TC}_TC.bam > ${result_path}${sample_name}_${TC}_sorted_TC.bam
  samtools index ${result_path}${sample_name}_${TC}_sorted_TC.bam
  bamCoverage -b ${result_path}${sample_name}_${TC}_sorted_TC.bam -o ${result_path}${sample_name}_${TC}_TC.bw

done




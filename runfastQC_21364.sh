#!/bin/bash
# QC with FASTQC
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
for ((lane=7;lane<=8;lane++))
do
for ((sample=1;sample<=16;sample++))
do
for ((pair=1;pair<=1;pair++))
do
bsub -q normal -G team152 -J pax.${lane}.${sample} -R "select[mem>1000] rusage[mem=1000]" -M1000 -o /lustre/scratch115/projects/paxgene/logs/21364_fastQC_${lane}_${sample}.out -e /lustre/scratch115/projects/paxgene/logs/21364_fastQC_${lane}.${sample}.err ${RNA_TOOLS}/FastQC/./fastqc -o /lustre/scratch115/projects/paxgene/fastQCres/ /lustre/scratch115/projects/paxgene/fastqFiles_21364/21364_${lane}.${sample}.${pair}.newcram.fastq
#/lustre/scratch115/projects/paxgene/fastqFiles_21364/21364_7.11.2.newcram.fastq
done
done
done




#!/bin/bash
# QC with FASTQC
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
for ((sample=1;sample<=16;sample++))
do
bsub -q normal -G team152 -J pax.${lane}.${sample} -R "select[mem>3000] rusage[mem=3000]" -M3000 -o /lustre/scratch115/projects/paxgene/logs/fastQC_${lane}_${sample}.out -e /lustre/scratch115/projects/paxgene/logs/fastQC_${lane}.${sample}.err ${RNA_TOOLS}/FastQC/./fastqc -o /lustre/scratch115/projects/paxgene/fastQCres/pastSortMeRNA/ /lustre/scratch115/projects/paxgene/fastqFiles/pastSortMeRNA/21121_${sample}.merged.fastq_non_rRNA.fastq
done

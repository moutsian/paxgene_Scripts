#!/bin/bash
# QC with FASTQC
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
for ((lane=7;lane<=8;lane++))
do
for ((sample=1;sample<=16;sample++))
do
for ((pair=2;pair<=2;pair++))
do
bsub -q normal -G team152 -J pax.${lane}.${sample} -R "select[mem>1000] rusage[mem=1000]" -M1000 -o /lustre/scratch115/projects/paxgene/logs/21121_fastQC_${lane}_${sample}.out -e /lustre/scratch115/projects/paxgene/logs/21121_fastQC_${lane}.${sample}.err ${RNA_TOOLS}/FastQC/./fastqc -o /lustre/scratch115/projects/paxgene/fastQCres/ /lustre/scratch115/projects/paxgene/fastqFiles/21121_${lane}_${sample}.${pair}.fastq
done
done
done




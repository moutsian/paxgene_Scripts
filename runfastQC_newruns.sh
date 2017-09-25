#!/bin/bash
# QC with FASTQC
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
run=$1
for ((lane=$2;lane<=$3;lane++))
do
for ((sample=$4;sample<=$5;sample++))
do
FILE="/lustre/scratch115/projects/paxgene/fastqfiles_newruns/${run}_${lane}_${sample}.1.sorted.fastq"
if [ ! -e "${FILE}" ];then
FILE="/lustre/scratch115/projects/paxgene/fastqfiles_newruns/${run}_${lane}.${sample}.1.sorted.fastq"
#echo ${FILE} changed 
fi
bsub -q normal -G team152 -J pax.${lane}.${sample}.${run} -R "select[mem>1000] rusage[mem=1000]" -M1000 -o /lustre/scratch115/projects/paxgene/logs/${run}_fastQC_${lane}_${sample}.out -e /lustre/scratch115/projects/paxgene/logs/${run}_fastQC_${lane}.${sample}.err ${RNA_TOOLS}/FastQC/./fastqc -o /lustre/scratch115/projects/paxgene/fastQCres_newruns/sorted/  ${FILE}
done
done




#!/bin/bash
# QC with FASTQC, after I unmerged the files again (post sortmerna)
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
for ((sample=1;sample<=16;sample++))
do
bsub -q normal -G team152 -J pax.fwd.${sample} -R "select[mem>3000] rusage[mem=3000]" -M3000 -o /lustre/scratch115/projects/paxgene/logs/fastQC_${sample}.fwd.out -e /lustre/scratch115/projects/paxgene/logs/fastQC.${sample}.fwd.err ${RNA_TOOLS}/FastQC/./fastqc -o /lustre/scratch115/projects/paxgene/fastQCres/pastTrimmomatic/ /lustre/scratch115/projects/paxgene/fastqFiles/pastTrimmomatic/21121.merged.${sample}.fwd.paired.postSortMeRNAtrimmo30.fastq
bsub -q normal -G team152 -J pax.fwd.${sample} -R "select[mem>3000] rusage[mem=3000]" -M3000 -o /lustre/scratch115/projects/paxgene/logs/fastQC_${sample}.bwd.out -e /lustre/scratch115/projects/paxgene/logs/fastQC.${sample}.bwd.err ${RNA_TOOLS}/FastQC/./fastqc -o /lustre/scratch115/projects/paxgene/fastQCres/pastTrimmomatic/ /lustre/scratch115/projects/paxgene/fastqFiles/pastTrimmomatic/21121.merged.${sample}.bwd.paired.postSortMeRNAtrimmo30.fastq
done

#!/bin/bash

#run Salmon
LIBRA=$1 #e.g. IU or ISR (which corresponds to the fr_firststrand option of TopHat
run=$2 #e.g. 22219, 22226, 22294
mkdir /lustre/scratch115/projects/paxgene/SalmonOutput_newruns/${LIBRA}
for ((lane=$3;lane<=$4;lane++));do
for ((sample=$5;sample<=$6;sample++));do
bsub -q normal -G team152 -J salmon.multi.${sample} -n4 -R "span[hosts=1] select[mem>7500] rusage[mem=7500]" -M7500 -o /lustre/scratch115/projects/paxgene/logs/${run}.${sample}.${LIBRA}.out -e /lustre/scratch115/projects/paxgene/logs/${run}.${sample}.${LIBRA}.err /software/team152/Salmon-0.8.2_linux_x86_64/bin/./salmon quant -i  /lustre/scratch115/projects/paxgene/new_index_files/GrCh38_cdna_from_ensembl/transcripts_index_0.8.2 -l "${LIBRA}" -1 /lustre/scratch115/projects/paxgene/fastqfiles_newruns/${run}_${lane}.${sample}.1.sorted.fastq  -2 /lustre/scratch115/projects/paxgene/fastqfiles_newruns/${run}_${lane}.${sample}.2.sorted.fastq  -o /lustre/scratch115/projects/paxgene/SalmonOutput_newruns/${LIBRA}/${run}_${lane}_${sample}_${LIBRA}.cdna_idx.sorted.raw -p 4
done
done


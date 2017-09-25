#!/bin/bash

# run MIXCR alignment on RNAseq data. Because here I am including partial assembling, this script 
# should be followed by assemblePartial. A new script is to be prepared for that.
# See https://mixcr.readthedocs.io/en/latest/rnaseq.html
LIBRA="21121"
for ((lane=7;lane<=8;lane++)); do
for ((sample=1;sample<=16;sample++));do
bsub -q normal -G team152 -J mixcr.${lane}.${sample} -n4 -R "span[hosts=1] select[mem>7500] rusage[mem=7500]" -M7500 -o /lustre/scratch115/projects/paxgene/logs/mixcr.${lane}.${sample}.out -e /lustre/scratch115/projects/paxgene/logs/mixcr.${lane}.${sample}.err /software/team152/mixcr-2.0.2/mixcr align -s hsa -t 4 --parameters rna-seq -OallowPartialAlignments=true -r /lustre/scratch115/projects/paxgene/logs/${lane}.${sample}.mixcr_align_report.txt /lustre/scratch115/projects/paxgene/fastqFiles/${LIBRA}_${lane}_${sample}.1.fastq /lustre/scratch115/projects/paxgene/fastqFiles/${LIBRA}_${lane}_${sample}.2.fastq /lustre/scratch115/projects/paxgene/tcrseq/${LIBRA}_${lane}_${sample}.alignments.vdjca
done ; done



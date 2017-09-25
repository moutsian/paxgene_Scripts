#!/bin/bash

# run MIXCR assemble on RNAseq data (step to be run after align and assemblePartial).
# See https://mixcr.readthedocs.io/en/latest/rnaseq.html
LIBRA="21121"
for ((lane=7;lane<=8;lane++)); do
for ((sample=1;sample<=16;sample++));do
bsub -q normal -G team152 -J mixcr.${lane}.${sample} -n4 -R "span[hosts=1] select[mem>7500] rusage[mem=7500]" -M7500 -o /lustre/scratch115/projects/paxgene/logs/mixcr.assemble.${lane}.${sample}.out -e /lustre/scratch115/projects/paxgene/logs/mixcr.assemble.${lane}.${sample}.err /software/team152/mixcr-2.0.2/mixcr assemble -r /lustre/scratch115/projects/paxgene/logs/${lane}.${sample}.mixcr_assemble_report.txt /lustre/scratch115/projects/paxgene/tcrseq/${LIBRA}_${lane}_${sample}.alignmentsRescued.vdjca /lustre/scratch115/projects/paxgene/tcrseq/${LIBRA}_${lane}_${sample}.clones.clns
done ; done



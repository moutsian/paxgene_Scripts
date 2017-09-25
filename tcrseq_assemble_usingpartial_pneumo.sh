#!/bin/bash

# run MIXCR assemble on RNAseq data (step to be run after align and assemblePartial).
# See https://mixcr.readthedocs.io/en/latest/rnaseq.html
LIBRA="19828"
for ((j=1;j<=7;j++)); do
for ((sample=1;sample<=8;sample++));do
bsub -q normal -G team152 -J mixcr.${j}.${sample} -n4 -R "span[hosts=1] select[mem>7500] rusage[mem=7500]" -M7500 -o /lustre/scratch115/projects/paxgene/logs/mixcr.assemble.pneumo${j}.${sample}.out -e /lustre/scratch115/projects/paxgene/logs/mixcr.assemble.pneumo${j}.${sample}.err /software/team152/mixcr-2.0.2/mixcr assemble -r /lustre/scratch115/projects/paxgene/logs/pneumo${j}.${sample}.mixcr_assemble_report.txt /lustre/scratch115/projects/paxgene/tcrseq/pneumo/${LIBRA}_${sample}_${j}.alignmentsRescued.vdjca /lustre/scratch115/projects/paxgene/tcrseq/pneumo/${LIBRA}_${sample}_${j}.clones.clns
done ; done



#!/bin/bash

# run MIXCR assemblePartial on RNAseq data.
# See https://mixcr.readthedocs.io/en/latest/rnaseq.html
#Also note that a second iteration is recommended in the MIXCR doumentation. Atm this is NOT being done in this script.
LIBRA="21121"
for ((lane=7;lane<=8;lane++)); do
for ((sample=1;sample<=16;sample++));do
bsub -q normal -G team152 -J mixcr.${lane}.${sample} -n4 -R "span[hosts=1] select[mem>7500] rusage[mem=7500]" -M7500 -o /lustre/scratch115/projects/paxgene/logs/mixcr.assemblePartial.${lane}.${sample}.out -e /lustre/scratch115/projects/paxgene/logs/mixcr.assemblePartial.${lane}.${sample}.err /software/team152/mixcr-2.0.2/mixcr assemblePartial -r /lustre/scratch115/projects/paxgene/logs/${lane}.${sample}.mixcr_assemblepartial_report.txt /lustre/scratch115/projects/paxgene/tcrseq/${LIBRA}_${lane}_${sample}.alignments.vdjca /lustre/scratch115/projects/paxgene/tcrseq/${LIBRA}_${lane}_${sample}.alignmentsRescued.vdjca
done ; done



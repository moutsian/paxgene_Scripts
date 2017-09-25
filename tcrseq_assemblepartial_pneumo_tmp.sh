#!/bin/bash

# run MIXCR assemblePartial on RNAseq data.
# See https://mixcr.readthedocs.io/en/latest/rnaseq.html
#Also note that a second iteration is recommended in the MIXCR doumentation. Atm this is NOT being done in this script.
LIBRA="19743"
for ((sample=1;sample<=1;sample++));do
for((j=3;j<=3;j++));do
bsub -q normal -G team152 -J mixcr.${j}.${sample} -n4 -R "span[hosts=1] select[mem>7500] rusage[mem=7500]" -M7500 -o /lustre/scratch115/projects/paxgene/logs/mixcr.assemblePartial.pneumo${j}.${sample}.out -e /lustre/scratch115/projects/paxgene/logs/mixcr.assemblePartial.pneumo${j}.${sample}.err /software/team152/mixcr-2.0.2/mixcr assemblePartial -r /lustre/scratch115/projects/paxgene/logs/pneumo${j}.${sample}.mixcr_assemblepartial_report.txt /lustre/scratch115/projects/paxgene/tcrseq/pneumo/${LIBRA}_${sample}_${j}.alignments.vdjca /lustre/scratch115/projects/paxgene/tcrseq/pneumo/${LIBRA}_${sample}_${j}.alignmentsRescued.vdjca
done ; done



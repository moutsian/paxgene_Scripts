#!/bin/bash

# run MIXCR assemble on RNAseq data (step to be run after align and assemblePartial).
# See https://mixcr.readthedocs.io/en/latest/rnaseq.html
LIBRA=$1
for ((j=1;j<=7;j++)); do
for ((sample=1;sample<=8;sample++));do
/software/team152/mixcr-2.0.2/mixcr exportClones -count -vGene -dGene -jGene -cGene -vFamily -dFamily -jFamily -cFamily -vHitScore -dHitScore -jHitScore -cHitScore /lustre/scratch115/projects/paxgene/tcrseq/pneumo/${LIBRA}_${sample}_${j}.clones.clns /lustre/scratch115/projects/paxgene/tcrseq/pneumo/${LIBRA}_${sample}_${j}.clones.txt
done ; done



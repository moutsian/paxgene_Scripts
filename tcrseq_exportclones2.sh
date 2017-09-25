#!/bin/bash

# run MIXCR assemble on RNAseq data (step to be run after align and assemblePartial).
# See https://mixcr.readthedocs.io/en/latest/rnaseq.html
LIBRA="21121"
for ((lane=7;lane<=8;lane++)); do
for ((sample=1;sample<=16;sample++));do
/software/team152/mixcr-2.0.2/mixcr exportClones -count -vGene -dGene -jGene -cGene -vFamily -dFamily -jFamily -cFamily -vHitScore -dHitScore -jHitScore -cHitScore /lustre/scratch115/projects/paxgene/tcrseq/${LIBRA}_${lane}_${sample}.clones.without_partial.clns /lustre/scratch115/projects/paxgene/tcrseq/${LIBRA}_${lane}_${sample}.clones.without_partial.txt
done ; done



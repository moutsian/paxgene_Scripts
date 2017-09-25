#!/bin/bash

# run MIXCR alignment on RNAseq data. Because here I am including partial assembling, this script
# should be followed by assemblePartial. A new script is to be prepared for that.
# See https://mixcr.readthedocs.io/en/latest/rnaseq.html


for entry in /lustre/scratch115/projects/pneumo-inf/19743_19828.fastqFiles.new.raw/*pair1*fastq
do
entry2=$(echo $entry|awk '{split($0,arr,"pair");{print arr[1]"pair2.fastq";}}')
#name=$(echo $entry|awk '{split($0,arr1,"raw");split(arr1[2],arr,"pair");{print arr[1];}}'|sed 's/#/_/g')
name=$(echo $entry|awk '{split($0,arr1,"raw/");split(arr1[2],arr,"pair");{print arr[1];}}'|sed 's/#/_/g')
echo ${name}
bsub -q normal -G team152 -J ${name} -n4 -R "span[hosts=1] select[mem>7500] rusage[mem=7500]" -M7500 -o /lustre/scratch115/projects/paxgene/tcrseq/pneumo/mixcr.${name}out -e /lustre/scratch115/projects/paxgene/tcrseq/pneumo/mixcr.${name}err /software/team152/mixcr-2.0.2/mixcr align -s hsa -t 4 --parameters rna-seq -OallowPartialAlignments=true -r /lustre/scratch115/projects/paxgene/tcrseq/pneumo/mixcr.${name}align.report ${entry} ${entry2} /lustre/scratch115/projects/paxgene/tcrseq/pneumo/${name}alignments.vdjca
done


#!/bin/bash

entry=$1;
entry2=$(echo $entry|awk '{split($0,arr,"pair");{print arr[1]"pair2.fastq";}}')
#name=$(echo $entry|awk '{split($0,arr,"pair");{print arr[1];}}'|sed 's/#/_/g')
name=$(echo $entry|awk '{split($0,arr1,"raw/");split(arr1[2],arr,"pair");{print arr[1];}}'|sed 's/#/_/g')
echo ${name}
/software/team152/mixcr-2.0.2/mixcr align -s hsa -t 4 --parameters rna-seq -OallowPartialAlignments=true -r /lustre/scratch115/projects/paxgene/tcrseq/pneumo/mixcr.${name}.align.report ${entry} ${entry2} /lustre/scratch115/projects/paxgene/tcrseq/pneumo/${name}alignments.vdjca


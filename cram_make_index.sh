#!/bin/bash
# SORTING BAM FILES
run=$1
for ((sample=1;sample<=16;sample++))
do
for ((lane=7;lane<=8;lane++))
do
/software/hgi/pkglocal/samtools-1.3.1-htslib-1.3.2/bin/samtools index /lustre/scratch115/projects/paxgene/sorted_bams_${run}/${run}_${lane}.${sample}.sorted.cram 
done
done

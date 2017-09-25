#!/bin/bash
# SORTING BAM FILES
run="21121"
for ((sample=1;sample<=16;sample++))
do
for ((lane=7;lane<=8;lane++))
do
bsub -q normal -G team152 -J ${run}${lane}.${sample} -R "select[mem>5000] rusage[mem=5000]" -M5000 -o /lustre/scratch115/projects/paxgene/logs/${run}_${lane}_${sample}.aligned.out -e /lustre/scratch115/projects/paxgene/logs/${run}_${lane}_${sample}.aligned.err /software/hgi/pkglocal/samtools-1.3.1-htslib-1.3.2/bin/samtools sort -o sorted_bams_${run}/${run}_${lane}.${sample}.sorted.cram /lustre/scratch115/projects/paxgene/cramfiles_${run}/${run}_${lane}.${sample}.cram
done
done

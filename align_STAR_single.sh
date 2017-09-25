#!/bin/bash
# take a fq_file's prefix as argument, and then, run STAR, there's some redundance here, caused modified from previous sh
#fq_list=$1
run=$1
lane=$2
sample=$3
n_threads=$4
#fq_indi_file=$fq_list
#fq_indi_prefix=`echo ${fq_indi_file} | awk -F/ '{print $NF}'`
fq_indi_prefix=${run}_${lane}_${sample}
out_dir=/lustre/scratch115/projects/paxgene/STAR_output/alignment/
whatever='.'
output_prefix=$out_dir$fq_indi_prefix$whatever
fq_dir=/lustre/scratch115/projects/paxgene/fastqfiles_newruns/
end1_suffix='.1.sorted.fastq'
end2_suffix='.2.sorted.fastq'
end1_fq=$fq_dir$fq_indi_prefix$end1_suffix
end2_fq=$fq_dir$fq_indi_prefix$end2_suffix

echo $end1_fq
echo $end2_fq

/software/team152/STAR-2.5.3a/bin/Linux_x86_64/STAR \
    --genomeDir /software/team152/ref/star_genome_index \
    --runThreadN $n_threads \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $output_prefix \
    --readFilesIn $end1_fq $end2_fq \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000


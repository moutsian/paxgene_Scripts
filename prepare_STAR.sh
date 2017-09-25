#!/bin/bash
# simple script, just to prepare some files before running STAR
ref=/lustre/scratch115/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa
anno=/software/team152/ref/gencode.v19.annotation.gtf.modified

/software/team152/STAR-2.5.3a/bin/Linux_x86_64/STAR \
    --runThreadN 20 \
    --runMode genomeGenerate \
    --genomeDir /software/team152/ref/star_genome_index/ \
    --genomeFastaFiles $ref \
    --sjdbGTFfile $anno \
    --sjdbOverhang 74

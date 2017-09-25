# BAM to FASTQ with samtools
run="21364"
for ((sample=1;sample<=16;sample++))
do
for ((lane=7;lane<=8;lane++))
do
bsub -q normal -G team152 -J $run$lane.$sample -R "select[mem>1000] rusage[mem=1000]" -M1000 -o /lustre/scratch115/projects/paxgene/logs/${run}_${lane}_${sample}.newcram.out -e /lustre/scratch115/projects/paxgene/logs/${run}_${lane}_${sample}.newcram.err /software/hgi/pkglocal/samtools-1.3.1-htslib-1.3.2/bin/samtools fastq -1 /lustre/scratch115/projects/paxgene/fastqFiles_${run}/${run}_${lane}_${sample}.1.newcram.fastq -2 /lustre/scratch115/projects/paxgene/fastqFiles_${run}/${run}_${lane}.${sample}.2.newcram.fastq /lustre/scratch115/projects/paxgene/cramfiles_${run}/${run}_${lane}.${sample}.cram
done
done

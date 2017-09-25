# BAM to FASTQ with samtools
run=$1
for ((sample=$4;sample<=$5;sample++))
do
for ((lane=$2;lane<=$3;lane++))
do
mkdir /lustre/scratch115/projects/paxgene/cramfiles_smartseq2/tmp${run}_${lane}_${sample}
bsub -q normal -G team152 -J $run$lane.$sample -R "select[mem>20000] rusage[mem=20000]" -M20000 -o /lustre/scratch115/projects/paxgene/logs/${run}_${lane}_${sample}.smartseq2.out -e /lustre/scratch115/projects/paxgene/logs/${run}_${lane}_${sample}.smartseq2.err /software/hgi/pkglocal/samtools-1.3.1-htslib-1.3.2/bin/samtools sort -m 8G -n  -o /lustre/scratch115/projects/paxgene/cramfiles_smartseq2/${run}_${lane}'#'${sample}.sorted.cram -T  /lustre/scratch115/projects/paxgene/cramfiles_smartseq2/tmp${run}_${lane}_${sample}/${run}_${lane}_${sample} /lustre/scratch115/projects/paxgene/cramfiles_smartseq2/${run}_${lane}'#'${sample}.cram 
#echo "submitted run ${run} for sample ${sample} in lanes ${lane}"
done
done

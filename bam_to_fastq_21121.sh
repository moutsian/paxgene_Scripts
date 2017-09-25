# BAM to FASTQ with samtools
run=$1
for ((sample=$4;sample<=$5;sample++))
do
for ((lane=$2;lane<=$3;lane++))
do
bsub -q normal -G team152 -J $run$lane.$sample -R "select[mem>1000] rusage[mem=1000]" -M1000 -o /lustre/scratch115/projects/paxgene/logs/${run}_${lane}_${sample}.${run}.out -e  /lustre/scratch115/projects/paxgene/logs/${run}_${lane}_${sample}.${run}.err /software/hgi/pkglocal/samtools-1.3.1-htslib-1.3.2/bin/samtools fastq -1  /lustre/scratch115/projects/paxgene/fastqfiles_newruns/${run}_${lane}.${sample}.1.sorted.fastq -2 /lustre/scratch115/projects/paxgene/fastqfiles_newruns/${run}_${lane}.${sample}.2.sorted.fastq /lustre/scratch115/projects/paxgene/cramfiles_${run}_fr_firststrand/${run}_${lane}'#'${sample}.sorted.cram
#echo "submitted run ${run} for sample ${sample} in lanes ${lane}"
done
done

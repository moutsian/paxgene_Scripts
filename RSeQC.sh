#align using STAR
run=$1
for ((lane=$2;lane<=$3;lane++));do
for ((sample=$4;sample<=$5;sample++));do
bsub -q normal -G team152 -J star.${sample} -n2 -R "span[hosts=1] select[mem>6000] rusage[mem=6000]" -M6000 -o /lustre/scratch115/projects/paxgene/logs/${run}.${sample}.RSeQC.out -e  /lustre/scratch115/projects/paxgene/logs/${run}.${sample}.RSeQC.err python /software/team152/RSeQC-2.6.4/scripts/read_duplication.py -i /lustre/scratch115/projects/paxgene/STAR_output/alignment/${run}_${lane}_${sample}.Aligned.sortedByCoord.out.bam -o  /lustre/scratch115/projects/paxgene/STAR_output/alignment/rseqc_stats/${run}_${lane}_${sample}
done
done


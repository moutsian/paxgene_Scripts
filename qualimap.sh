#align using STAR
run=$1
for ((lane=$2;lane<=$3;lane++));do
for ((sample=$4;sample<=$5;sample++));do
# a) cluster version - this will produce the statistics but not the html report and the plots
#bsub -q normal -G team152 -J star.${sample} -n6 -R "span[hosts=1] select[mem>8000] rusage[mem=8000]" -M8000 -o /lustre/scratch115/projects/paxgene/logs/${run}.${sample}.qualimap.out -e  /lustre/scratch115/projects/paxgene/logs/${run}.${sample}.qualimap.err /software/team152/qualimap_v2.2.1/./qualimap bamqc -bam /lustre/scratch115/projects/paxgene/STAR_output/alignment/${run}_${lane}_${sample}.Aligned.sortedByCoord.out.bam --java-mem-size=6G -nt 6
#b) standard version
/software/team152/qualimap_v2.2.1/./qualimap bamqc -bam /lustre/scratch115/projects/paxgene/STAR_output/alignment/${run}_${lane}_${sample}.Aligned.sortedByCoord.out.bam --java-mem-size=6G -nt 6

done
done


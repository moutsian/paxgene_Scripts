#align using STAR
run=$1
for ((lane=$2;lane<=$3;lane++));do
for ((sample=$4;sample<=$5;sample++));do
bsub -q normal -G team152 -J star.${sample} -n6 -R "span[hosts=1] select[mem>34000] rusage[mem=34000]" -M34000 -o /lustre/scratch115/projects/paxgene/logs/${run}.${sample}.star.out -e /lustre/scratch115/projects/paxgene/logs/${run}.${sample}.star.err sh align_STAR_single.sh "$run" "$lane" "$sample" 6
done
done


#run SortmeRNA on merged fastq files - if these are past sortmeRNA, ensure you have first unpaired them using the provided script by the sortMERNA ppl.
for ((sample=1;sample<=16;sample++));do
# for ((sample=3;sample<=3;sample++));do
bsub -q normal -G team152 -J trimmo.${sample} -n4 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -M10000 -o /lustre/scratch115/projects/paxgene/logs/trimmo.${sample}.out -e /lustre/scratch115/projects/paxgene/logs/trimmo.${sample}.err java -jar /software/team152/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 4 /lustre/scratch115/projects/paxgene/fastqFiles/pastSortMeRNA/21121_${sample}.merged.fastq_non_rRNA.1.fastq /lustre/scratch115/projects/paxgene/fastqFiles/pastSortMeRNA/21121_${sample}.merged.fastq_non_rRNA.2.fastq /lustre/scratch115/projects/paxgene/fastqFiles/pastTrimmomatic/21121.merged.${sample}.fwd.paired.postSortMeRNAtrimmo30.fastq.gz /lustre/scratch115/projects/paxgene/fastqFiles/pastTrimmomatic/21121.merged.${sample}.fwd.unpaired.postSortMeRNAtrimmo30.fastq.gz /lustre/scratch115/projects/paxgene/fastqFiles/pastTrimmomatic/21121.merged.${sample}.bwd.paired.postSortMeRNAtrimmo30.fastq.gz /lustre/scratch115/projects/paxgene/fastqFiles/pastTrimmomatic/21121.merged.${sample}.bwd.unpaired.postSortMeRNAtrimmo30.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30
done



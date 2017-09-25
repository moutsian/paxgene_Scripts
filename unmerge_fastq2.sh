#merge fastq files to run SortmeRNA - after we've merged both lanes into a single fastq
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
for ((sample=1;sample<=16;sample++));do
bsub -q normal -G team152 -J unmerge.${sample} -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -M5000 -o /lustre/scratch115/projects/paxgene/logs/21121.unmerge_${sample}.out -e /lustre/scratch115/projects/paxgene/logs/21121.unmerge.${sample}.err  sh /software/team152/sortmerna-2.1-linux-64/scripts/unmerge-paired-reads.sh /lustre/scratch115/projects/paxgene/fastqFiles/pastSortMeRNA/21121_${sample}.merged.fastq_non_rRNA.fastq /lustre/scratch115/projects/paxgene/fastqFiles/pastSortMeRNA/21121_${sample}.merged.fastq_non_rRNA.1.fastq /lustre/scratch115/projects/paxgene/fastqFiles/pastSortMeRNA/21121_${sample}.merged.fastq_non_rRNA.2.fastq
done



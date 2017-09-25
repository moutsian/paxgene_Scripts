#merge fastq files to run SortmeRNA
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
for ((lane=7;lane<=8;lane++)); do
for ((sample=1;sample<=16;sample++));do
bsub -q normal -G team152 -J merge.${lane}.${sample} -n1 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o /lustre/scratch115/projects/paxgene/logs/21121.merge${lane}_${sample}.out -e /lustre/scratch115/projects/paxgene/logs/21121.merge.${lane}.${sample}.err sh /software/team152/sortmerna-2.1-linux-64/scripts/merge-paired-reads.sh /lustre/scratch115/projects/paxgene/fastqFiles/21121_${lane}_${sample}.1.fastq /lustre/scratch115/projects/paxgene/fastqFiles/21121_${lane}_${sample}.2.fastq /lustre/scratch115/projects/paxgene/fastqFiles/21121_${lane}_${sample}.merged.fastq 
done ; done



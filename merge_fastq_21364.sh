#merge fastq files to run SortmeRNA
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
dataset="21364"
for ((lane=7;lane<=8;lane++)); do
for ((sample=1;sample<=16;sample++));do
bsub -q normal -G team152 -J merge.${lane}.${sample} -n1 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o /lustre/scratch115/projects/paxgene/logs/${dataset}.merge.${lane}_${sample}.out -e  /lustre/scratch115/projects/paxgene/logs/${dataset}.merge.${lane}.${sample}.err sh /software/team152/sortmerna-2.1-linux-64/scripts/merge-paired-reads.sh /lustre/scratch115/projects/paxgene/fastqFiles_${dataset}/${dataset}_${lane}.${sample}.1.newcram.fastq /lustre/scratch115/projects/paxgene/fastqFiles_${dataset}/${dataset}_${lane}.${sample}.2.newcram.fastq /lustre/scratch115/projects/paxgene/fastqFiles_${dataset}/${dataset}_${lane}_${sample}.merged.fastq 
#/lustre/scratch115/projects/paxgene/fastqFiles_21364/21364_8.8.1.newcram.fastq
done ; done



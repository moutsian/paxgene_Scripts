#merge fastq files to run SortmeRNA
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
dataset=$1 # e.g. 21364

for ((lane=$2;lane<=$3;lane++)); do
for ((sample=$4;sample<=$5;sample++));do
bsub -q normal -G team152 -J merge.${lane}.${sample} -n1 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o /lustre/scratch115/projects/paxgene/logs/${dataset}.merge.${lane}_${sample}.out -e /lustre/scratch115/projects/paxgene/logs/${dataset}.merge.${lane}.${sample}.err sh /software/team152/sortmerna-2.1-linux-64/scripts/merge-paired-reads.sh /lustre/scratch115/projects/paxgene/fastqfiles_newruns/${dataset}_${lane}_${sample}.1.fastq /lustre/scratch115/projects/paxgene/fastqfiles_newruns/${dataset}_${lane}_${sample}.2.fastq /lustre/scratch115/projects/paxgene/fastqfiles_newruns/${dataset}_${lane}_${sample}.merged.fastq 
#/lustre/scratch115/projects/paxgene/fastqFiles_21364/21364_8.8.1.newcram.fastq
done ; done



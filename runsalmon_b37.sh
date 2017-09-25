#run Salmon
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
for ((lane=7;lane<=8;lane++)); do
for ((sample=1;sample<=16;sample++));do
bsub -q normal -G team152 -J salmon.${lane}.${sample} -n4 -R "span[hosts=1] select[mem>7500] rusage[mem=7500]" -M7500 -o /lustre/scratch115/projects/paxgene/SalmonOutput/21121.b37.${lane}.${sample}.out -e /lustre/scratch115/projects/paxgene/logs/21121.b37.${lane}.${sample}.err ${RNA_TOOLS}Salmon-0.7.2_linux_x86_64/bin/./salmon quant -i /lustre/scratch113/teams/anderson/users/jga/UsefulFiles/Homo_sapiens.GRCh37.rel75.cdna.all.Salmon.idx -l A -1 /lustre/scratch115/projects/paxgene/fastqFiles/21121_${lane}_${sample}.1.fastq -2 /lustre/scratch115/projects/paxgene/fastqFiles/21121_${lane}_${sample}.2.fastq -o /lustre/scratch115/projects/paxgene/SalmonOutput/21121_${lane}.${sample}.b37.raw -p 4
done ; done



#run Salmon
LIBRA=$1
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
for ((sample=1;sample<=16;sample++));do
bsub -q normal -G team152 -J salmon.multi.${sample} -n4 -R "span[hosts=1] select[mem>7500] rusage[mem=7500]" -M7500 -o /lustre/scratch115/projects/paxgene/logs/21121.${sample}.${LIBRA}.multi.out -e /lustre/scratch115/projects/paxgene/logs/21121.${sample}.${LIBRA}.multi.err ${RNA_TOOLS}Salmon-0.7.2_linux_x86_64/bin/./salmon quant -i /lustre/scratch115/projects/paxgene/new_index_files/GRCh38_15_plus_hs38d1.known.salmon.idx -l "${LIBRA}" -1 /lustre/scratch115/projects/paxgene/fastqFiles/21121_7_${sample}.1.fastq /lustre/scratch115/projects/paxgene/fastqFiles/21121_8_${sample}.1.fastq -2 /lustre/scratch115/projects/paxgene/fastqFiles/21121_7_${sample}.2.fastq /lustre/scratch115/projects/paxgene/fastqFiles/21121_8_${sample}.2.fastq -o /lustre/scratch115/projects/paxgene/SalmonOutput/21121.${sample}.new_index.${LIBRA}.multi.raw -p 4
done



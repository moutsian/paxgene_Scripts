#run Salmon
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
LIBRA="IU"
for ((sample=1;sample<=16;sample++));do
bsub -q normal -G team152 -J salmon.ni.${sample} -n4 -R "span[hosts=1] select[mem>7500] rusage[mem=7500]" -M7500 -o /lustre/scratch115/projects/paxgene/logs/21121.${sample}.ni.${LIBRA}.out -e /lustre/scratch115/projects/paxgene/logs/21121.${sample}.ni.${LIBRA}.err ${RNA_TOOLS}Salmon-0.7.2_linux_x86_64/bin/./salmon quant -i /lustre/scratch115/projects/paxgene/new_index_files/GRCh38_15_plus_hs38d1.known.salmon.idx -l A -1 /lustre/scratch115/projects/paxgene/fastqFiles/21121_${sample}.1.fastq -2 /lustre/scratch115/projects/paxgene/fastqFiles/21121_${sample}.2.fastq -o /lustre/scratch115/projects/paxgene/SalmonOutput/21121.${sample}.new_index.mergedbyme.${LIBRA}.raw -p 4
done



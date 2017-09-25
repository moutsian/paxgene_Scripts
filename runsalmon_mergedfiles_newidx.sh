#run Salmon with the new idx I built using the same transcriptome ref as the one mentioned at the bamfile
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
for ((sample=1;sample<=16;sample++));do
bsub -q normal -G team152 -J salmon.${sample} -n4 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -M10000 -o /lustre/scratch115/projects/paxgene/logs/21121.${sample}.newidx.isr.out -e /lustre/scratch115/projects/paxgene/logs/21121.${sample}.isr.err ${RNA_TOOLS}Salmon-0.7.2_linux_x86_64/bin/./salmon quant -i /lustre/scratch115/projects/paxgene/new_index_files/GRCh38_15_plus_hs38d1.known.salmon.idx -l ISR -1 /lustre/scratch115/projects/paxgene/fastqFiles/pastTrimmomatic/21121.merged.${sample}.fwd.paired.postSortMeRNAtrimmo30.fastq -2 /lustre/scratch115/projects/paxgene/fastqFiles/pastTrimmomatic/21121.merged.${sample}.bwd.paired.postSortMeRNAtrimmo30.fastq -o /lustre/scratch115/projects/paxgene/SalmonOutput_after_filtering/21121.${sample}.postfilter.newidx.ISR.raw -p 4
done



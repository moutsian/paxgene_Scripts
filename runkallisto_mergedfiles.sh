#run Kallisto
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
for ((sample=1;sample<=16;sample++));do
bsub -q normal -G team152 -J kallisto.${sample} -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -M10000 -o /lustre/scratch115/projects/paxgene/logs/21121.${sample}.kallisto.out -e /lustre/scratch115/projects/paxgene/logs/21121.${sample}.kallisto.err ${RNA_TOOLS}/kallisto_linux-v0.43.0/kallisto quant -i /lustre/scratch113/teams/anderson/users/jga/UsefulFiles/Homo_sapiens.GRCh38.rel79.cdna.all.Kallisto.idx -o /lustre/scratch115/projects/paxgene/KallistoOutput_after_filtering/21121.${sample}.postfilter.raw -b 100 /lustre/scratch115/projects/paxgene/fastqFiles/pastTrimmomatic/21121.merged.${sample}.fwd.paired.postSortMeRNAtrimmo30.fastq -2 /lustre/scratch115/projects/paxgene/fastqFiles/pastTrimmomatic/21121.merged.${sample}.bwd.paired.postSortMeRNAtrimmo30.fastq -o 
done



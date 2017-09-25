#run Salmon
RNA_TOOLS="/lustre/scratch113/teams/anderson/users/jga/ToolBox/RNA_tools/"
for ((lane=7;lane<=8;lane++)); do
for ((sample=1;sample<=16;sample++));do
fragratio=$(grep compatible_fragment_ratio /lustre/scratch115/projects/paxgene/SalmonOutput/21121_${lane}.${sample}.raw/lib_format_counts.json)
echo ${lane} ${sample} ${fragratio}
done ; done



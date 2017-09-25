#merge the samples from each lane into one.
#merge across lanes:
for ((sample=1;sample<=16;sample++));do
echo "sample: ${sample}"
 cat /lustre/scratch115/projects/paxgene/fastqFiles/21121_7_${sample}.1.fastq /lustre/scratch115/projects/paxgene/fastqFiles/21121_8_${sample}.1.fastq > /lustre/scratch115/projects/paxgene/fastqFiles/21121_${sample}.1.fastq
cat /lustre/scratch115/projects/paxgene/fastqFiles/21121_7_${sample}.2.fastq /lustre/scratch115/projects/paxgene/fastqFiles/21121_8_${sample}.2.fastq >  /lustre/scratch115/projects/paxgene/fastqFiles/21121_${sample}.2.fastq
done





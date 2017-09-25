#!/bin/bash
#for FILE in SalmonOutput_newruns/IU_sorted/*sorted.raw/ ; do
#for FILE in /lustre/scratch115/projects/paxgene/SalmonOutput_newruns/IU/22891_*cdna_idx.sorted.raw/;do
for FILE in /lustre/scratch115/projects/paxgene/SalmonOutput_newruns/ISR_21121_sorted/*/;do
echo $FILE
paste /lustre/scratch115/projects/paxgene/SalmonOutput/IU_21364/21364.3.new_index.IU.multi.raw/quant.sf  ${FILE}/quant.sf|awk '{print $1,$7,$8,$9,$10}' > ${FILE}/quant.sf.updated
cp ${FILE}/quant.sf  ${FILE}/quant.sf.original
cat ${FILE}/quant.sf |awk '{split($1,arr,".");print arr[1],$2,$3,$4,$5}' >  ${FILE}/quant.sf.updated
mv ${FILE}/quant.sf.updated ${FILE}/quant.sf
sed -i s'/ /\t/g' ${FILE}/quant.sf
done

#!/bin/bash
for FILE in SalmonOutput_newruns/IU_sorted/*sorted.raw/ ; do
echo $FILE
paste SalmonOutput/IU_21364/21364.3.new_index.IU.multi.raw/quant.sf  ${FILE}/quant.sf|awk '{print $1,$7,$8,$9,$10}' > ${FILE}/quant.sf.updated
#mv ${FILE}/quant.sf.updated ${FILE}/quant.sf
#sed -i s'/ /\t/g' ${FILE}/quant.sf
done

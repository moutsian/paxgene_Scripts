#!/bin/bash
libra=$1;
for((sample=1;sample<=8;i++))
do
for((lane=7;lane<=8;lane++));do
sh align_STAR_single.sh "$libra" "$lane" "$sample"
done;done

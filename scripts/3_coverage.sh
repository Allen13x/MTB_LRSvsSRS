#!/bin/bash
mkdir Coverage
cd Coverage
for j in $(ls -1 ../reads)
do
for i in $(cat ../sample_list)
do
mosdepth -t 8 -x -b ../REF/mtb_region1k.bed -T 1,2,4,6,8,12,16,20,40 ${j}_${i}_1k ../reads/${j}/Bam/${i}_*.bam
mosdepth -t 8 -x -b ../REF/pe_ppe.bed -T 1,2,4,6,8,12,16,20,40 ${j}_${i}_pe-ppe ../reads/${j}/Bam/${i}_*.bam
done
done 

#!/bin/bash

ls -1 reads/LRS/ | cut -f1 -d '_' > sample_list
for sample in $(cat sample_list)
do

Ratatosk -s reads/SRS/${sample}_*R1* -s reads/SRS/${sample}_*R2* -l reads/LRS/${sample}_*.fastq.gz -o reads/Hybrid/${sample}_Hybrid.fastq.gz

done

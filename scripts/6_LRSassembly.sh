#!/bin/bash

mkdir Assembly
mkdir Assembly/LRS
mkdir Assembly/Hybrid

for i in $(ls -1 reads/LRS/)
do
sample=$(echo $i | sed 's/_.*//g')
flye --nano-raw reads/LRS/${i} --out-dir Assembly/LRS/${sample}  -i 3 --genome-size 4.4m

cp Assembly/LRS/${sample}/assembly.fasta Assembly/LRS/LRS_${sample}.fasta

done


for i in $(ls -1 reads/Hybrid/)
do
sample=$(echo $i | sed 's/_.*//g')
flye --nano-raw reads/Hybrid/${i} --out-dir Assembly/Hybrid/${sample}  -i 3 --genome-size 4.4m

cp Assembly/Hybrid/${sample}/assembly.fasta Assembly/Hybrid/Hybrid_${sample}.fasta

done



#!/bin/bash

mkdir Assembly
mkdir Assembly/LRS
mkdir Assembly/Hybrid

# LRS Assembly

for i in $(ls -1 reads/LRS/)
do
sample=$(echo $i | sed 's/_.*//g')
flye --nano-raw reads/LRS/${i} --out-dir Assembly/LRS/${sample}  -i 3 --genome-size 4.4m

cp Assembly/LRS/${sample}/assembly.fasta Assembly/LRS/LRS_${sample}.fasta

done

# Hybrid Reads Assembly
for i in $(ls -1 reads/Hybrid/)
do
sample=$(echo $i | sed 's/_.*//g')
flye --nano-raw reads/Hybrid/${i} --out-dir Assembly/Hybrid/${sample}  -i 3 --genome-size 4.4m

cp Assembly/Hybrid/${sample}/assembly.fasta Assembly/Hybrid/Hybrid_${sample}.fasta

done



mkdir Assembly
mkdir Assembly/SRS
mkdir Assembly/HybA

# SRS and Hybrid Assembly

for sample in $(ls -1 reads/SRS/ | sed 's/_.*//g' | sort | uniq)
do
unicycler -1 reads/SRS/${sample}_*R1* -2 reads/SRS/${sample}_*R2* -o Assembly/SRS/${sample}
cp Assembly/SRS/${sample}/assembly.fasta Assembly/SRS/SRS_${sample}.fasta


unicycler -1 reads/SRS/${sample}_*R1* -2 reads/SRS/${sample}_*R2* -l reads/LRS/${sample}_* -o Assembly/HybA/${sample}

cp Assembly/HybA/${sample}/assembly.fasta Assembly/HybA/HybA_${sample}.fasta

done

# Quast analysis

quast -r REF/M._tuberculosis_H37Rv_2015-11-13.fasta -g REF/M._tuberculosis_H37Rv_2015-11-13_genes.bed -o Assembly/quast Assembly/LRS/LRS_*.fasta Assembly/Hybrid/Hybrid_*.fasta Assembly/HybA/HybA_*.fasta Assembly/SRS/SRS_*.fasta






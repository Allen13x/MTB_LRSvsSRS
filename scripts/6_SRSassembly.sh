#!/bin/bash

mkdir Assembly
mkdir Assembly/SRS
mkdir Assembly/HybA

for sample in $(ls -1 reads/SRS/ | sed 's/_.*//g' | sort | uniq)
do
unicycler -1 reads/SRS/${sample}_*R1* -2 reads/SRS/${sample}_*R2* -o Assembly/SRS/${sample}
cp Assembly/SRS/${sample}/assembly.fasta Assembly/SRS/SRS_${sample}.fasta


unicycler -1 reads/SRS/${sample}_*R1* -2 reads/SRS/${sample}_*R2* -l reads/LRS/${sample}_* -o Assembly/HybA/${sample}

cp Assembly/HybA/${sample}/assembly.fasta Assembly/HybA/HybA_${sample}.fasta

done


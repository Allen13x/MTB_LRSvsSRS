#!/bin/bash


Ratatosk -s reads/SRS/${sample}_*R1* -s reads/SRS/${sample}_*R2* -l reads/LRS/${sample}*.fastq.gz -o reads/Hybrid/${sample}

# Long and short reads sequencing: two inclusive technologies for the hybrid investigation of the Mycobacterium tuberculosis genome

![image](https://user-images.githubusercontent.com/72440375/186920399-67b26432-1154-41c9-8f1d-12c6b86ada78.png)


## Settings

1) Use the scripts/0_req.txt to build the conda environment for the analysis


2) Create a folder reads with three subfolders: LRS, SRS and Hybrid and put your .fastq files in the respective folders 
3) Install MTBseq from  https://github.com/ngs-fzb/MTBseq_source
4) Install Ratatosk from https://github.com/DecodeGenetics/Ratatosk
5) Install R with the packages tidyverse and rstatix



## Generate Hybrid reads

From the main directory launch the script 1_Hybrid_generate-reads.sh

## Running MTBseq

From the main directory launch the three scripts:
2_Hyb_mtbseq.sh
2_LRS_mtbseq.sh
2_SRS_mtbseq.sh

## Genome Coverage

From the main directory launch the script: 3a_coverage.sh

Open R and and use the code in 3b_Coverage.R to compare the breadth coverage at 8x with the different approaches

## Variant calling and Comparative Analysis

The variants were already called with the MTBseq step. They can be found in the relative folder (e.g reads/SRS/Called).
From the main directory launch the script: 4_tree.sh
This will produce a folder Trees containing the transmission trees of the samples according the different approaches.

For visualization use the site: https://achtman-lab.github.io/GrapeTree/MSTree_holder.html

## Drug resistance

For each of the approaches enter the folder reads/*/Called.
From there launch the script 5_WHO-resistance.py

The results will be in the file qcr26_all_full.csv

## De Novo Assembly

From the main irectory launch the script 6a_Assembly.sh

The assemblies will be in the folder Assembly

For comparison and visualization use the R script 6b_assembly.R

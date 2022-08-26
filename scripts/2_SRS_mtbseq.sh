#!/bin/bash 




for i in *_R2*.fastq.gz; do basename=`ls $i | cut -d "_" -f 1`;x=$(date +%d%m%y) ; echo $x ; mv $i ${basename}_SRS-Q4-RP1-PH0_150bp_R2.fastq.gz; done
for i in *_R1*.fastq.gz; do basename=`ls $i | cut -d "_" -f 1`;x=$(date +%d%m%y) ; echo $x ; mv $i ${basename}_SRS-Q4-RP1-PH0_150bp_150bp_R1.fastq.gz; done
for i in *R1.fastq.gz; do x=` ls $i | cut -d "_" -f1,2` ; echo $x >> samplelist ; done
perl -p -e "s/_/\t/g" samplelist > sample_list



MTBseq --step TBbwa --mincovf 1 --mincovr 1 --lowfreq_vars --minfreq 5 --minphred20 1 --samples sample_list --threads 8
MTBseq --step TBrefine --mincovf 1 --mincovr 1 --lowfreq_vars --minfreq 5 --minphred20 1 --samples sample_list --threads 8
MTBseq --step TBpile --mincovf 1 --mincovr 1 --lowfreq_vars --minfreq 5 --minphred20 1 --samples sample_list --threads 8
MTBseq --step TBlist --mincovf 1 --mincovr 1 --lowfreq_vars --minfreq 5 --minphred20 1 --samples sample_list --threads 8
MTBseq --step TBvariants --mincovf 1 --mincovr 1 --lowfreq_vars --minfreq 5 --minphred20 1 --samples sample_list --threads 8
MTBseq --step TBvariants --samples sample_list --threads 8 --minbqual 4 --categories ../REF/MTB_Gene_Categories.txt --minphred 0
MTBseq --step TBstats --samples sample_list --threads 8  --minbqual 4 --categories ../REF/MTB_Gene_Categories.txt --minphred 0
MTBseq --step TBstrains --sample sample_list --threads 8  --minbqual 4 --categories ../REF/MTB_Gene_Categories.txt --minphred 0


A=`date +"%d%b%y%a%H%M%S"`

MTBseq --step TBjoin --continue --samples sample_list --project proj${A}  --minbqual 4 --categories ../REF/MTB_Gene_Categories.txt --minphred 0

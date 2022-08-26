#!/bin/bash 

cd reads/LRS


for i in *.fastq.gz; do basename=`ls $i | cut -d "." -f 1 | cut -d "_" -f 1`; echo $x; mv $i ${basename}_LRS-Q4-RP1-PH0_R1.fastq.gz; done
for i in *.fastq.gz; do x=` ls $i | cut -d "_" -f1,2` ; echo $x >> samplelist; done
cat samplelist | tr "_" "\t" > sample_list

mkdir Bam
mkdir GATK_Bam
mkdir Mpileup
mkdir Position_Tables
mkdir Called
mkdir Statistics
mkdir Classification

for sample in $(cat sample_list)
do
bwa mem -t 8 /beegfs/datasets/buffer/ric.cirillo/MTB/M._tuberculosis_H37Rv_2015-11-13.fasta ${sample}*.fastq.gz > Bam/${sample}.sam 2>> Bam/${sample}.bamlog
samtools view -@ 8 -b -T /beegfs/datasets/buffer/ric.cirillo/MTB/M._tuberculosis_H37Rv_2015-11-13.fasta -o Bam/${sample}.bam Bam/${sample}.sam 2>> Bam/${sample}.bamlog
samtools sort -@ 8 -T /tmp/${sample}.sorted -o Bam/${sample}.sorted.bam Bam/${sample}.bam 2>> Bam/${sample}.bamlog
samtools index -b Bam/${sample}.sorted.bam 2>> Bam/${sample}.bamlog
samtools rmdup Bam/${sample}.sorted.bam Bam/${sample}.nodup.bam 2>> Bam/${sample}.bamlog
samtools index -b Bam/${sample}.nodup.bam 2>> Bam/${sample}.bamlog



rm Bam/${sample}.sam Bam/${sample}.bam Bam/${sample}.sorted.bam Bam/${sample}.sorted.bam.bai

mv Bam/${sample}.nodup.bam Bam/${sample}_nBP.bam
mv Bam/${sample}.nodup.bam.bai Bam/${sample}_nBP.bam.bai


mkdir temp_Bam

cat <(samtools view -H Bam/${sample}_nBP.bam) <(paste <(samtools view Bam/${sample}_nBP.bam | cut -f1-10 ) <(samtools view Bam/${sample}_nBP.bam | cut -f 11 | tr "$(cat ../../REF/ascii_string)" "K")) | samtools view -b -o temp_Bam/${sample}_dump.bam -
picard AddOrReplaceReadGroups I=temp_Bam/${sample}_dump.bam O=temp_Bam/${sample}_final.bam RGPU=unit1 RGID=11 RGLB=LaneX RGSM=AnySampleName RGPL=illumina


samtools index temp_Bam/${sample}_final.bam
gatk3 -Xmx30g --analysis_type RealignerTargetCreator --reference_sequence ../../REF/M._tuberculosis_H37Rv_2015-11-13.fasta --input_file temp_Bam/${sample}_final.bam --downsample_to_coverage 10000 --num_threads 8 --out GATK_Bam/${sample}.gatk.intervals 2>> GATK_Bam/${sample}.gatk.bamlog

gatk3 -Xmx30g --analysis_type IndelRealigner --reference_sequence ../../REF/M._tuberculosis_H37Rv_2015-11-13.fasta --input_file temp_Bam/${sample}_final.bam --defaultBaseQualities 4 --targetIntervals GATK_Bam/${sample}.gatk.intervals --noOriginalAlignmentTags --out GATK_Bam/${sample}.realigned.bam 2>> GATK_Bam/${sample}.gatk.bamlog

gatk3 -Xmx30g --analysis_type BaseRecalibrator --reference_sequence ../../REF/M._tuberculosis_H37Rv_2015-11-13.fasta --input_file GATK_Bam/${sample}.realigned.bam --knownSites ../../REF/MTB_Base_Calibration_List.vcf --maximum_cycle_value 400000  --num_cpu_threads_per_data_thread 8 --out GATK_Bam/${sample}.gatk.grp 2>>GATK_Bam/${sample}.gatk.bamlog

gatk3 -Xmx30g -T --analysis_type PrintReads --reference_sequence ../../REF/M._tuberculosis_H37Rv_2015-11-13.fasta --input_file GATK_Bam/${sample}.realigned.bam --BQSR GATK_Bam/${sample}.gatk.grp --num_cpu_threads_per_data_thread 8 --out GATK_Bam/${sample}_nBP.gatk.bam 2>> GATK_Bam/${sample}.gatk.bamlog


samtools index GATK_Bam/${sample}_nBP.gatk.bam

rm GATK_Bam/*.realigned.*

samtools mpileup -B -A -x -Q ${minbqual} -f /opt/common/tools/ric.cosr/miniconda3/envs/mtbseq/share/mtbseq-1.0.4-1/var/ref/M._tuberculosis_H37Rv_2015-11-13.fasta -o Mpileup/${sample}_nBP.gatk.mpileup GATK_Bam/${sample}_nBP.gatk.bam


rm -r temp_Bam




MTBseq --step TBlist --mincovf 4 --mincovr 4 --minfreq 75  --minbqual 4 --categories ../../REF/MTB_Gene_Categories.txt --minphred 0 --samples sample_list --threads 8
MTBseq --step TBvariants --mincovf 1 --mincovr 1 --lowfreq_vars --minfreq 5 --minphred20 1 --categories ../../REF/MTB_Gene_Categories.txt
MTBseq --step TBvariants  --mincovf 4 --mincovr 4 --minfreq 75  --minbqual 4 --categories ../../REF/MTB_Gene_Categories.txt --minphred 0 --samples sample_list --threads 8
#MTBseq --step TBvariants --samples sample_list --threads 8
MTBseq --step TBstats  --minbqual 4 --mincovf 4 --mincovr 4 --minfreq 75 --categories ../../REF/MTB_Gene_Categories.txt --minphred 0 --samples sample_list --threads 8
MTBseq --step TBstrains  --minbqual 4 --mincovf 4 --mincovr 4 --minfreq 75 --categories ../../REF/MTB_Gene_Categories.txt --minphred 0 --sample sample_list --threads 8


A=`date +"%d%b%y%a%H%M%S"`

MTBseq --step TBjoin --continue --samples sample_list --project proj${A}  --minbqual 4 --categories ../../REF/MTB_Gene_Categories.txt --minphred 0





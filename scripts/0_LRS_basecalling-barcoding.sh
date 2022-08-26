#!/bin/bash

guppy_basecaller --input_path fast5/ --save_path fastq_barcode_hac/ --disable_pings --compress_fastq --progress_stats_frequency 60 --config dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 10 -x auto --barcode_kits SQK-RBK004

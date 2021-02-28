#!/bin/bash

. ~/.profile

echo "running on $HOSTNAME"
echo "sample name = $1"

start=$(date +%s)
echo "start time: $start"


in_dir=/path/to/COVIRT_all_data_rRNA_removed/02-AlignedData
gatk_path=/path/to/hsp/COVID19/anaconda3/envs/COVIRT_GATK/bin
ref=/path/to/RNAseq_Metrics_refs


echo "gatk version: "
$gatk_path/gatk --version

echo "SAMPLE: $1"

call='$gatk_path/gatk CollectRnaSeqMetrics -I $in_dir/$1/$1_Aligned.sortedByCoord_sorted.out.bam -O $in_dir/$1/$1_Metrics --REF_FLAT $ref/Homo_sapiens.GRCh38_and_Sars_cov_2.ASM985889v3_refFlat.txt -STRAND SECOND_READ_TRANSCRIPTION_STRAND --METRIC_ACCUMULATION_LEVEL ALL_READS --RIBOSOMAL_INTERVALS $ref/Homo_sapiens.GRCh38.100_and_SARS-CoV-2_rRNA_interval_list.txt'

echo $call
eval $call

end=$(date +%s)
echo "end time: $end"
runtime_s=$(echo $(( end - start )))
echo "total run time(s): $runtime_s"
sec_per_min=60
sec_per_hr=3600
runtime_m=$(echo "scale=2; $runtime_s / $sec_per_min;" | bc)
echo "total run time(m): $runtime_m"
runtime_h=$(echo "scale=2; $runtime_s / $sec_per_hr;" | bc)
echo "total run time(h): $runtime_h"


#!/bin/bash

. ~/.profile

echo "running on $HOSTNAME"
echo "sample name = $1"

start=$(date +%s)
echo "start time: $start"

in_dir=/path/to/00-RawData/Fastq
trimlog_dir=/path/to/01-TrimmedData/Trimming_Reports
out_dir=/path/to/01-TrimmedData/Fastq
adapter_dir=/path/to/trimmomatic/adapters

call="trimmomatic PE \
-threads 4 \
-phred33 \
-trimlog $trimlog_dir/$1_trimming.log \
-summary $trimlog_dir/$1_trimming_summary.txt \
-validatePairs \
$in_dir/$1.R1.fastq.gz \
$in_dir/$1.R2.fastq.gz \
$out_dir/$1_R1_P_trimmed.fq.gz \
$out_dir/$1_R1_U_trimmed.fq.gz \
$out_dir/$1_R2_P_trimmed.fq.gz \
$out_dir/$1_R2_U_trimmed.fq.gz \
ILLUMINACLIP:$adapter_dir/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
LEADING:20 \
TRAILING:20 \
MINLEN:15"

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

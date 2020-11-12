#!/bin/bash

. ~/.profile

echo "running on $HOSTNAME"
echo "sample name = $1"

start=$(date +%s)
echo "start time: $start"

in_dir=/path/to/01-TrimmedData/Fastq
out_dir=/path/to/01-TrimmedData/Fastq_RG

mkdir $out_dir/$1

module load /path/to/hsp/COVID19/anaconda3.modulefile
cd /path/to/hsp/COVID19/gdc-fastq-splitter
source ./path/to/bin/activate

echo "gdc-fastq-splitter version:"
gdc-fastq-splitter --version

call="gdc-fastq-splitter -o $out_dir/$1/$1_ $in_dir/$1_R1_P_trimmed.fq.gz $in_dir/$1_R2_P_trimmed.fq.gz"

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

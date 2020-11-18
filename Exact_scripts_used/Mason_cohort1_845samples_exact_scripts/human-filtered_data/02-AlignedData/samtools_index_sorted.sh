#!/bin/bash

. ~/.profile

echo "running on $HOSTNAME"
echo "sample name = $1"

start=$(date +%s)
echo "start time: $start"

in_dir=/path/to/02-AlignedData_RG

echo "samtools version: "
samtools --version

call="samtools index -@ 4 $in_dir/$1/$1_Aligned.sortedByCoord_sorted.out.bam"

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

#!/bin/bash

. ~/.profile

echo "running on $HOSTNAME"
echo "sample name = $2"
echo "chr = chr$1"

start=$(date +%s)
echo "start time: $start"

gatk_path=/path/to/COVID19/anaconda3/envs/COVIRT_GATK/bin
in_dir=/path/to/02-AlignedData
inter_dir=/path/to/GATK_analysis/03-GATKIntermediates
out_dir=/path/to/GATK_analysis/04-GATKOutput
genome_ref=/path/to/genome_files/Homo_sapiens/ensembl_release100
vcf_ref=/path/to/genome_files/Homo_sapiens/ensembl_release100/vcf_files

mkdir $out_dir/$2/chr$1


call6='$gatk_path/gatk --java-options "-Xmx100G" HaplotypeCaller -R $genome_ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I $inter_dir/$2/$2_BSQR-applied.out.bam -O $out_dir/$2/chr$1/$2_chr$1.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20 -G StandardAnnotation -G AS_StandardAnnotation --intervals $1'

echo $call6
eval $call6


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

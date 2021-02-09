#!/usr/bin/env bash

split_reads_dir="../Fastq_Input_Files_for_Testing/split_reads"

sample="ENVHA_332764270_HA_all-reads"

# TODO: Remove basename when in workdir
flowcells=`ls -1 ${split_reads_dir}/split_${sample}_*_*_R*.fq.gz \
    | xargs -L 1 basename \
    | sed s/split_// \
    | sed s/${sample}// \
    | cut -d "_" -f 2 \
    | uniq`

for flowcell in ${flowcells}; do
    echo ${flowcell}
    lanes=`ls -1 ${split_reads_dir}/split_${sample}_${flowcell}_*_R*.fq.gz \
    	      | xargs -L 1 basename \
	      | sed s/split_// \
	      | sed s/${sample}// \
	      | cut -d "_" -f 3 \
	      | sort \
	      | uniq`
    echo ${lanes}

    outSAMattrRGlineStr=""
    readOneFilesIn=""
    for lane in ${lanes}; do

	str="ID:${flowcell}.${lane} PL:ILLUMINA PU:${flowcell}.${lane} LB:${sample} SM:${sample}"
	if [ -z "${outSAMattrRGlineStr}" ]; then
	    outSAMattrRGlineStr="$str"
	else
	    outSAMattrRGlineStr="${outSAMattrRGlineStr}, $str"
	fi
	
	str=split_${sample}_${flowcell}_${lane}_R1.fq.gz
	if [ -z "${readOneFilesIn}" ]; then
	    readOneFilesIn="${str}"
	else
	    readOneFilesIn="${readOneFilesIn},${str}"
	fi

    done
    readTwoFilesIn=`echo ${readOneFilesIn} | sed s/R1/R2/g`
    echo "${outSAMattrRGlineStr}"
    echo "${readOneFilesIn}"
    echo "${readTwoFilesIn}"
    echo

done

# for flowcell in ${flowcells}; do
#     echo ${flowcell}
# 
#     reads_one_files=`ls -1 ${split_reads_dir}/split_${sample}_${flowcell}_*_R1.fq.gz \
#         | xargs -L 1 basename \
#         | tr "\n" ","`
#     echo ${reads_one_files}
# 
#     reads_two_files=`echo ${reads_one_files} | sed s/R1/R2/g`
#     echo ${reads_two_files}
# 
#     echo
# 
# done

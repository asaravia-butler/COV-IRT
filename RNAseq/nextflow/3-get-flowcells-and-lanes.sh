#!/usr/bin/env bash

split_reads_dir="../Fastq_Input_Files_for_Testing/split_reads"

sample="ENVHA_332764270_HA_all-reads"

# TODO: Remove basename when in workdir
flowcells=`ls -1 ${split_reads_dir}/split_${sample}_* \
    | xargs -L 1 basename \
    | sed s/split_// \
    | sed s/${sample}// \
    | cut -d "_" -f 2 \
    | uniq`

echo ${flowcells}

for flowcell in ${flowcells}; do
    echo ${flowcell}

    reads_one_files=`ls -1 ${split_reads_dir}/split_${sample}_${flowcell}_*_R1.fq.gz \
        | xargs -L 1 basename \
        | tr "\n" ","`
    echo ${reads_one_files}

    reads_two_files=`echo ${reads_one_files} | sed s/R1/R2/g`
    echo ${reads_two_files}

    echo

done

#!/usr/bin/env nextflow

// TODO: Understand data storage on NASA cluster

params.raw_reads_dir = "/home/ubuntu/COV-IRT/RNAseq/Fastq_Input_Files_for_Testing"
params.ensembl_data_dir = "/home/ubuntu/COV-IRT-Data"

process createSampleList {

  publishDir params.raw_reads_dir, mode: "copy"

  output:
    file "samples.txt" into samples_file_ch

  """
  ls ${params.raw_reads_dir}/*.R1.fastq.gz \
    | xargs basename \
    | sed s/.R1.fastq.gz// \
    > samples.txt
  """
}

raw_reads_files_ch = Channel.fromPath(params.raw_reads_dir + "/*.fastq.gz")

params.raw_fastqc_dir = params.raw_reads_dir + "/raw_fastqc"

process createRawReadsQC {

  label "align"

  publishDir params.raw_fastqc_dir, mode: "copy"

  input:
    // TODO: Use where needed
    // file 'samples.txt' from samples_file_ch
    file raw_reads_file from raw_reads_files_ch

  output:
    file "*" into raw_fastqc_ch

  """
  fastqc ${raw_reads_file}
  """
}

params.raw_multiqc_dir = params.raw_reads_dir + "/raw_multiqc"

process compileRawReadsQC {

  label "align"

  publishDir params.raw_multiqc_dir, mode: "copy"

  input:
    file "*" from raw_fastqc_ch.collect()

  output:
    file "*" into raw_multiqc_ch

  """
  multiqc -n raw_multiqc -f ${params.raw_fastqc_dir}
  """
}

raw_reads_file_pairs_ch = Channel.fromFilePairs(params.raw_reads_dir + "/*.R{1,2}.fastq.gz")

// TODO: Understand setting to use system maximum threads
params.numberOfThreads = 2
params.trimmed_reads_dir = params.raw_reads_dir + "/trimmed_reads"

process trimRawReads {

  label "align"
  
  publishDir params.trimmed_reads_dir, mode: "copy"

  input:
    set sample, file(raw_reads_file_pair) from raw_reads_file_pairs_ch

  output:
    file "*_P_trimmed.fq.gz" into trimmed_reads_files_ch
    file "*_R1_P_trimmed.fq.gz" into trimmed_reads_one_files_ch
    file "*_R2_P_trimmed.fq.gz" into trimmed_reads_two_files_ch

  """
  trimmomatic PE \
    -threads ${params.numberOfThreads} \
    -phred33 \
    -trimlog ${sample}_trimming.log \
    -summary ${sample}_trimming_summary.txt \
    -validatePairs \
    ${raw_reads_file_pair} \
    ${sample}_R1_P_trimmed.fq.gz \
    ${sample}_R1_U_trimmed.fq.gz \
    ${sample}_R2_P_trimmed.fq.gz \
    ${sample}_R2_U_trimmed.fq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
    LEADING:20 \
    TRAILING:20 \
    MINLEN:15
  """
}

params.trimmed_fastqc_dir = params.trimmed_reads_dir + "/trimmed_fastqc"

process createTrimmedReadsQC {

  label "align"

  publishDir params.trimmed_fastqc_dir, mode: "copy"

  input:
    file trimmed_reads_file from trimmed_reads_files_ch

  output:
    file "*" into trimmed_fastqc_ch

  """
  fastqc ${trimmed_reads_file}
  """
}

params.trimmed_multiqc_dir = params.trimmed_reads_dir + "/trimmed_multiqc"

process compileTrimmedReadsQC {

  label "align"

  publishDir params.trimmed_multiqc_dir, mode: "copy"

  input:
    file "*" from trimmed_fastqc_ch.collect()

  output:
    file "*" into trimmed_multiqc_ch

  """
  multiqc -n trimmed_multiqc -f ${params.trimmed_fastqc_dir}
  """
}

process splitTrimmedReads {

  input:
    file trimmed_reads_one_file from trimmed_reads_one_files_ch
    file trimmed_reads_two_file from trimmed_reads_two_files_ch

  """
  sample_one=`echo ${trimmed_reads_one_file} | sed s/_R1_P_trimmed.fq.gz//`
  sample_two=`echo ${trimmed_reads_two_file} | sed s/_R2_P_trimmed.fq.gz//`
  # TODO: Decide if this test is required
  # Ensure sample reads file pairs are not interleaved
  if [ \${sample_one} != \${sample_two} ]; then
    exit 1
  fi
  docker run -v ${params.trimmed_reads_dir}:/opt --rm \
    quay.io/kmhernan/gdc-fastq-splitter -o split_\${sample_one}_ \
    ${trimmed_reads_one_file} ${trimmed_reads_two_file}
  """
}

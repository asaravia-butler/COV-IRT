#!/usr/bin/env nextflow

params.raw_reads_dir = "/Users/raymondleclair/Projects/IQT/COV-IRT/RNAseq/Fastq_Input_Files_for_Testing"

process listSamples {
  """
  cd ${params.raw_reads_dir}
  ls -1 *.R1.fastq.gz | sed s/.R1.fastq.gz// > samples.txt
  """
}

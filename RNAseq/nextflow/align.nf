#!/usr/bin/env nextflow

raw_reads_files_ch = Channel.fromPath(params.raw_reads_dir + "/*.fastq.gz")

params.raw_fastqc_dir = params.raw_reads_dir + "/raw_fastqc"

process createRawReadsQC {

  label "COVIRT_fastq_to_alignment"

  publishDir params.raw_fastqc_dir, mode: "copy"

  input:
    file raw_reads_file from raw_reads_files_ch

  output:
    file "*" into raw_fastqc_ch

  """
  fastqc ${raw_reads_file}
  """
}

params.raw_multiqc_dir = params.raw_reads_dir + "/raw_multiqc"

process compileRawReadsQC {

  label "COVIRT_fastq_to_alignment"

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

// TODO: Automate setting of this value
params.numberOfThreads = 16

params.trimmed_reads_dir = params.raw_reads_dir + "/trimmed_reads"

process trimRawReads {

  label "COVIRT_fastq_to_alignment"
  
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

  label "COVIRT_fastq_to_alignment"

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

  label "COVIRT_fastq_to_alignment"

  publishDir params.trimmed_multiqc_dir, mode: "copy"

  input:
    file "*" from trimmed_fastqc_ch.collect()

  output:
    file "*" into trimmed_multiqc_ch

  """
  multiqc -n trimmed_multiqc -f ${params.trimmed_fastqc_dir}
  """
}

params.split_reads_dir = params.raw_reads_dir + "/split_reads"

process splitTrimmedReads {

  // TODO: Specify Docker image?

  publishDir params.split_reads_dir, mode: "copy"

  input:
    file trimmed_reads_one_file from trimmed_reads_one_files_ch
    file trimmed_reads_two_file from trimmed_reads_two_files_ch

  output:
    env sample into split_reads_sample_ch
    file "split_*_R1.fq.gz" into split_reads_one_files_ch
    file "split_*_R2.fq.gz" into split_reads_two_files_ch
    file "split_*.json" into split_reads_json_files_ch

  """
  sample=`echo ${trimmed_reads_one_file} | sed s/_R1_P_trimmed.fq.gz//`
  # Remove links and copy in reads files
  rm ${trimmed_reads_one_file}
  rm ${trimmed_reads_two_file}
  cp ${params.trimmed_reads_dir}/${trimmed_reads_one_file} .
  cp ${params.trimmed_reads_dir}/${trimmed_reads_two_file} .

  # Split reads
  docker run -v \${PWD}:/opt --rm \
    quay.io/kmhernan/gdc-fastq-splitter -o split_\${sample}_ \
    ${trimmed_reads_one_file} ${trimmed_reads_two_file}
  """
}

params.filtered_reads_dir = params.raw_reads_dir + "/filtered_reads"

process filterSplitReads {

  label "COVIRT_HTStream"

  publishDir params.filtered_reads_dir, mode: "copy"

  input:
    env sample from split_reads_sample_ch
    file split_reads_one_file from split_reads_one_files_ch
    file split_reads_two_file from split_reads_two_files_ch

  output:
    env sample into filtered_reads_sample_ch
    file "split_*_R1.fastq.gz" into filtered_reads_one_files_ch
    file "split_*_R2.fastq.gz" into filtered_reads_two_files_ch
    file "split_*.log" into filtered_reads_log_files_ch

  """
  for reads_one_file in `ls split_\${sample}_*_R1.fq.gz`; do
    reads_two_file=`echo \${reads_one_file} | sed s/R1/R2/`
    base_name=`echo \${reads_one_file} \
      | sed s/_R1.fq.gz//`
    # TODO: Need to get the reference file to the worker
    hts_SeqScreener \
      -L \${base_name}_htsStats.log \
      -1 ${split_reads_one_file} \
      -2 ${split_reads_two_file} \
      -s ${params.COVIRT_home}/RNAseq/Reference_Files/Hsapiens_rRNA_RefSeq_seq_w_mitrRNA_ITS_ETS.fasta \
      -x 0.20 \
      -f \${base_name}

  done
  """
}

params.aligned_reads_dir = params.raw_reads_dir + "/aligned_reads"

process alignReadsToReferenceGenome {

  label "COVIRT_fastq_to_alignment"

  publishDir params.aligned_reads_dir, mode: "copy"

  input:
    env sample from filtered_reads_sample_ch
    file split_reads_json_file from split_reads_json_files_ch

  output:
    file "*ReadsPerGene.out.tab" into reads_counts_files_ch
    file "*_Aligned.sortedByCoord.out.bam" into aligned_reads_files_ch
    file "*_Aligned.toTranscriptome.out.bam" into aligned_transcriptome_reads_files_ch

  """
  # Copy in reads files (can't use links anyway)
  cp ${params.filtered_reads_dir}/split_\${sample}_*.fastq.gz .

  # Handle multiple flowcells per sample
  outSAMattrRGline=""
  readOneFilesIn=""
  readTwoFilesIn=""
  for flowcell in `ls -1 split_\${sample}_*_R1.fastq.gz \
    | xargs -L 1 basename \
    | sed s/split_// \
    | sed s/\${sample}// \
    | cut -d "_" -f 2 \
    | uniq`; do

    # Handle multiple lanes per flowcell
    for lane in `ls -1 split_\${sample}_\${flowcell}_*_R1.fastq.gz \
      | xargs -L 1 basename \
      | sed s/split_// \
      | sed s/\${sample}// \
      | cut -d "_" -f 3 \
      | uniq`; do

      # Construct the read file paths
      if [ -n "\${readOneFilesIn}" ]; then
        readOneFilesIn="\${readOneFilesIn},"
      fi
      readOneFilesIn="\${readOneFilesIn}split_\${sample}_\${flowcell}_\${lane}_R1.fastq.gz"
      if [ -n "\${readTwoFilesIn}" ]; then
        readTwoFilesIn="\${readTwoFilesIn},"
      fi
      readTwoFilesIn="\${readTwoFilesIn}split_\${sample}_\${flowcell}_\${lane}_R2.fastq.gz"

      # Find the multiplex barcode
      barcode=`grep "multiplex_barcode" split_\${sample}_\${flowcell}_\${lane}_R1.report.json \
        | cut -d ":" -f 2 \
        | sed s/\\"//g \
        | sed s/,//g \
        | xargs`

      # Construct the read group line
      if [ -n "\${outSAMattrRGline}" ]; then
        outSAMattrRGline="\${outSAMattrRGline} , "
      fi
      outSAMattrRGline="\${outSAMattrRGline}ID:\${flowcell}.\${lane} PL:ILLUMINA PU:\${flowcell}.\${lane}.\${barcode} LB:\${sample} SM:\${sample}"
    done
  done
  readFilesIn="\${readOneFilesIn} \${readTwoFilesIn}"

  # Align reads
  STAR \
    --twopassMode Basic \
    --limitBAMsortRAM 65000000000 \
    --genomeDir ${params.genome_dir} \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS nM NM MD jM jI MC ch \
    --outSAMattrRGline \${outSAMattrRGline} \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --readFilesCommand zcat \
    --runThreadN ${params.numberOfThreads} \
    --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
    --chimOutJunctionFormat 1 \
    --chimSegmentMin 20 \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --outFileNamePrefix \${sample}_\${flowcell}_ \
    --readFilesIn \${readFilesIn}
  """
}

process generateStarCountsTable {

  publishDir params.aligned_reads_dir, mode: "copy"

  input:
    file "*" from reads_counts_files_ch.collect()

  output:
    file "*" into counts_table_file_ch

  """
  #!/usr/bin/env r

  print("Make STAR counts table")
  print("")

  ### Import data
  ff <- list.files(
    file.path("${params.aligned_reads_dir}"),
    pattern="ReadsPerGene.out.tab", full.names=TRUE
  )

  ## Remove the first 4 lines
  counts.files <- lapply(ff, read.table, skip = 4)

  ## Get counts aligned to the second, reverse, strand
  counts <- as.data.frame(sapply(counts.files, function(x) x[ , 4 ]))

  ## Get sample names
  samples <- sub(
    "_[a-zA-Z0-9]*_[0-9]*_ReadsPerGene.out.tab", "",
    sub("split_", "", basename(ff))
  )

  ## Add column and row names
  colnames(counts) <- samples
  row.names(counts) <- counts.files[[1]]\$V1

  ### Export unnormalized counts table
  write.csv(counts, file='STAR_Unnormalized_Counts.csv')

  print("Session Info below: ")
  print("")
  sessionInfo()
  """
}

// TODO: Automate setting of these values
params.AvailableMemoryPerThread = "1G"
params.NumberOfThreads = 16

process sortAndIndexAlignedReads {

  label "COVIRT_fastq_to_alignment"

  publishDir params.aligned_reads_dir, mode: "copy"

  input:
    file aligned_reads_file from aligned_reads_files_ch

  output:
    file "*" into sorted_reads_files_ch

  """
  sorted_reads_file=`echo ${aligned_reads_file} | sed s/.out.bam/_sorted.out.bam/`
  samtools sort \
    -m ${params.AvailableMemoryPerThread} --threads ${params.NumberOfThreads} \
    -o \${sorted_reads_file} ${aligned_reads_file}
  samtools index \
    -@ ${params.NumberOfThreads} \${sorted_reads_file}
  """
}

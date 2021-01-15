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
    | xargs basename -a \
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
params.numberOfThreads = 16
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

params.split_reads_dir = params.raw_reads_dir + "/split_reads"

process splitTrimmedReads {

  input:
    file trimmed_reads_one_file from trimmed_reads_one_files_ch
    file trimmed_reads_two_file from trimmed_reads_two_files_ch

  // TODO: Need and output block
  
  """
  sample_one=`echo ${trimmed_reads_one_file} | sed s/_R1_P_trimmed.fq.gz//`
  sample_two=`echo ${trimmed_reads_two_file} | sed s/_R2_P_trimmed.fq.gz//`
  # TODO: Decide if this test is required
  # Ensure sample reads file pairs are not interleaved
  if [ \${sample_one} != \${sample_two} ]; then
    exit 1
  fi
  mkdir -p ${params.split_reads_dir}
  # Docker doesn't like symbolic links
  cp ${params.trimmed_reads_dir}/${trimmed_reads_one_file} ${params.split_reads_dir}
  cp ${params.trimmed_reads_dir}/${trimmed_reads_two_file} ${params.split_reads_dir}
  docker run -v ${params.split_reads_dir}:/opt --rm \
    quay.io/kmhernan/gdc-fastq-splitter -o split_\${sample_one}_ \
    ${trimmed_reads_one_file} ${trimmed_reads_two_file}
  rm ${params.split_reads_dir}/${trimmed_reads_one_file}
  rm ${params.split_reads_dir}/${trimmed_reads_two_file}
  """
}

split_reads_file_pairs_ch = Channel.fromFilePairs(params.split_reads_dir + "/*_R{1,2}.fq.gz")

params.genome_dir = params.ensembl_data_dir + "/filtered"
params.aligned_reads_dir = params.raw_reads_dir + "/aligned_reads"
// TODO: What value should this be?
params.index = 3

process alignReadsToReferenceGenome {

  label "align"

  publishDir params.aligned_reads_dir, mode: "copy"

  input:
    set sample, file(split_reads_file_pair) from split_reads_file_pairs_ch

  output:
    file "*" into aligned_reads_ch

  """
  rm ${split_reads_file_pair[0]}
  rm ${split_reads_file_pair[1]}
  # TODO: Why is this needed?
  cp ${params.split_reads_dir}/${split_reads_file_pair[0]} .
  cp ${params.split_reads_dir}/${split_reads_file_pair[1]} .
  STAR \
    --twopassMode Basic \
    --limitBAMsortRAM 65000000000 \
    --genomeDir ${params.genome_dir} \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS nM NM MD jM jI MC ch \
    --outSAMattrRGline \
        ID:flowcell.laneX PL:ILLUMINA PU:flowcell.laneX.${params.index} LB:${sample} SM:${sample}, \
        ID:flowcell.laneY PL:ILLUMINA PU:flowcell.laneY.${params.index} LB:${sample} SM:${sample} \
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
    --outFileNamePrefix ${sample}_ \
    --readFilesIn ${split_reads_file_pair}
  """
}

process generateStarCountsTable {

  publishDir params.aligned_reads_dir, mode: "copy"

  input:
    file 'samples.txt' from samples_file_ch
    file "*" from aligned_reads_ch.collect()

  output:
    file "*" into counts_table_ch

  """
  #!/usr/bin/env r

  print("Make STAR counts table")
  print("")

  # TODO: Get sample names from ff
  ### Pull in sample names
  study <- read.csv(Sys.glob(file.path("${params.raw_reads_dir}", "samples.txt")),
    header=FALSE, row.names=1, stringsAsFactors=TRUE)

  ### Import data
  ff <- list.files(file.path("${params.aligned_reads_dir}"),
    pattern="ReadsPerGene.out.tab", full.names=TRUE)

  ## Remove the first 4 lines
  counts.files <- lapply(ff, read.table, skip = 4)

  ## Get counts aligned to the second, reverse, strand
  counts <- as.data.frame(sapply(counts.files, function(x) x[ , 4 ]))

  ## Add column and row names
  colnames(counts) <- rownames(study)
  row.names(counts) <- counts.files[[1]]\$V1

  ### Export unnormalized counts table
  write.csv(counts, file='STAR_Unnormalized_Counts.csv')

  ### print session info
  print("Session Info below: ")
  print("")
  sessionInfo()
  """
}

aligned_reads_files_ch = Channel.fromPath(params.aligned_reads_dir + "/*_Aligned.sortedByCoord.out.bam")

params.AvailableMemoryPerThread = "1G"
params.NumberOfThreads = 16

process sortAndIndexAlignedReads {

  label "align"

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

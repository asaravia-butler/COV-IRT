#!/usr/bin/env nextflow

raw_reads_files_ch = Channel.fromPath(params.raw_reads_dir + "/*.fastq.gz")

params.raw_fastqc_dir = params.COVIRT_Code + "/00-RawData/FastQC_Reports"

process createRawReadsQC {

    label "COVIRT_fastq_to_alignment"

    publishDir params.raw_fastqc_dir, mode: "copy"

    input:
    file raw_reads_file from raw_reads_files_ch

    output:
    file "*" into raw_fastqc_ch

    """
    fastqc -t 4 ${raw_reads_file}
    """
}

params.raw_multiqc_dir = params.raw_fastqc_dir + "/raw_multiqc_report"

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
params.numberOfThreads = 10

params.trimmed_data_dir = params.COVIRT_Code + "/01-TrimmedData"

params.trimmed_reads_dir = params.trimmed_data_dir + "/Fastq"

params.trimmed_logs_dir = params.trimmed_data_dir + "/Trimming_Reports"

process trimRawReads {

    label "COVIRT_fastq_to_alignment"
  
    publishDir params.trimmed_reads_dir, mode: "copy"

    input:
    set sample, file(raw_reads_file_pair) from raw_reads_file_pairs_ch

    output:
    file "*_P_trimmed.fq.gz" into trimmed_P_reads_files_ch
    file "*_U_trimmed.fq.gz" into trimmed_U_reads_files_ch
    file "*_R1_P_trimmed.fq.gz" into trimmed_reads_one_files_ch
    file "*_R2_P_trimmed.fq.gz" into trimmed_reads_two_files_ch

    """
    trimmomatic PE \
        -threads ${params.numberOfThreads} \
        -phred33 \
        -trimlog ${params.trimmed_logs_dir}/${sample}_trimming.log \
        -summary ${params.trimmed_logs_dir}/${sample}_trimming_summary.txt \
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

params.trimmed_fastqc_dir = params.trimmed_data_dir + "/FastQC_Reports"

params.trimmed_P_fastqc_dir = params.trimmed_fastqc_dir + "/paired_reads"


process createTrimmedPReadsQC {

    label "COVIRT_fastq_to_alignment"

    publishDir params.trimmed_P_fastqc_dir, mode: "copy"

    input:
    file trimmed_P_reads_file from trimmed_P_reads_files_ch

    output:
    file "*" into trimmed_P_fastqc_ch

    """
    fastqc -t 4 ${trimmed_P_reads_file}
    """
}


params.trimmed_U_fastqc_dir = params.trimmed_fastqc_dir + "/unpaired_reads"


process createTrimmedUReadsQC {

    label "COVIRT_fastq_to_alignment"

    publishDir params.trimmed_U_fastqc_dir, mode: "copy"

    input:
    file trimmed_U_reads_file from trimmed_U_reads_files_ch

    output:
    file "*" into trimmed_U_fastqc_ch

    """
    fastqc -t 4 ${trimmed_U_reads_file}
    """
}


params.trimmed_P_multiqc_dir = params.trimmed_fastqc_dir + "/trimmed_P_multiqc_report"

process compileTrimmedPReadsQC {

    label "COVIRT_fastq_to_alignment"

    publishDir params.trimmed_P_multiqc_dir, mode: "copy"

    input:
    file "*" from trimmed_P_fastqc_ch.collect()

    output:
    file "*" into trimmed_P_multiqc_ch

    """
    multiqc -n trimmed_P_multiqc -f ${params.trimmed_P_fastqc_dir}
    """
}


params.trimmed_U_multiqc_dir = params.trimmed_fastqc_dir + "/trimmed_U_multiqc_report"

process compileTrimmedUReadsQC {

    label "COVIRT_fastq_to_alignment"

    publishDir params.trimmed_U_multiqc_dir, mode: "copy"

    input:
    file "*" from trimmed_U_fastqc_ch.collect()

    output:
    file "*" into trimmed_U_multiqc_ch

    """
    multiqc -n trimmed_U_multiqc -f ${params.trimmed_U_fastqc_dir}
    """
}


params.split_reads_dir = params.trimmed_data_dir + "/Fastq_RG"

process splitTrimmedReads {

    label "gdc_fastq_splitter"

    publishDir params.split_reads_dir, mode: "copy"

    input:
    file trimmed_reads_one_file from trimmed_reads_one_files_ch
    file trimmed_reads_two_file from trimmed_reads_two_files_ch

    output:
    env sample into split_reads_sample_ch
    file "*/split_*_R1.fq.gz" into split_reads_one_files_ch
    file "*/split_*_R2.fq.gz" into split_reads_two_files_ch
    file "*/split_*.json" into split_reads_json_files_ch

  """
  sample=`echo ${trimmed_reads_one_file} | sed s/_R1_P_trimmed.fq.gz//`
  echo Sample \$sample

  mkdir \${sample}

  # Remove links and copy in reads files
  # rm ${trimmed_reads_one_file}
  # rm ${trimmed_reads_two_file}
  # cp ${params.trimmed_reads_dir}/${trimmed_reads_one_file} .
  # cp ${params.trimmed_reads_dir}/${trimmed_reads_two_file} .

  # Split reads
  # docker run -v \${PWD}:/opt --rm \
  # quay.io/kmhernan/
  gdc-fastq-splitter -o \${sample}/split_\${sample}_ \
    ${trimmed_reads_one_file} ${trimmed_reads_two_file}
  """
}

params.filtered_reads_dir = params.trimmed_data_dir + "/Fastq_RG_NOrRNA"

params.htstream_logs_dir = params.trimmed_data_dir + "/HTStream_logs"

process filterSplitReads {

    label "COVIRT_HTStream"

    publishDir params.filtered_reads_dir, mode: "copy", pattern: "*/*.fastq.gz"
    publishDir params.htstream_logs_dir, mode: "copy", pattern: "*_htsStats.log"

    input:
    env sample from split_reads_sample_ch
    file split_reads_one_file from split_reads_one_files_ch
    file split_reads_two_file from split_reads_two_files_ch

    output:
    env sample into filtered_reads_sample_ch
    file "*/split_*_R1.fastq.gz" into filtered_reads_one_files_ch
    file "*/split_*_R2.fastq.gz" into filtered_reads_two_files_ch
    file "*/split_*.log" into filtered_reads_log_files_ch

    """
    mkdir \${sample}

    # Handle multiple flowcells per sample
    for flowcell in `ls -1 \${sample}/split_\${sample}_*_R1.fq.gz \
        | xargs -L 1 basename \
        | sed s/split_// \
        | sed s/\${sample}// \
        | cut -d "_" -f 2 \
        | uniq`; do

        # Handle multiple lanes per flowcell
        for lane in `ls -1 \${sample}/split_\${sample}_\${flowcell}_*_R1.fq.gz \
            | xargs -L 1 basename \
            | sed s/split_// \
            | sed s/\${sample}// \
            | cut -d "_" -f 3 \
            | uniq`; do
            reads_one_file="\${sample}/split_\${sample}_\${flowcell}_\${lane}_R1.fq.gz"
            reads_two_file="\${sample}/split_\${sample}_\${flowcell}_\${lane}_R2.fq.gz"
            base_name="\${sample}/split_\${sample}_\${flowcell}_\${lane}"
            # TODO: Need to get the reference file to the worker
            hts_SeqScreener \
                -L \${base_name}_htsStats.log \
                -1 \${reads_one_file} \
                -2 \${reads_two_file} \
                -s ${params.COVIRT_Data}/rRNA_refs/Hsapiens_rRNA_RefSeq_sequence.fasta \
                -x 0.20 \
                -f \${sample}/\${base_name}
        done
    done
    """
}

params.aligned_reads_dir = params.COVIRT_Code + "/02-AlignedData"

process alignReadsToReferenceGenome {

    label "COVIRT_fastq_to_alignment"

    publishDir params.aligned_reads_dir, mode: "copy"

    input:
    env sample from filtered_reads_sample_ch
    file split_reads_json_file from split_reads_json_files_ch

    output:
    file "*/*ReadsPerGene.out.tab" into reads_counts_files_ch
    file "*/*_Aligned.sortedByCoord.out.bam" into aligned_reads_files_ch
    file "*/*_Aligned.toTranscriptome.out.bam" into aligned_transcriptome_reads_files_ch
    file "*/*_Log.final.out" into star_log_files_ch

    """
    mkdir \${sample}

    # Copy in reads files (can't use links anyway)
    cp ${params.filtered_reads_dir}/\${sample}/split_\${sample}_*.fastq.gz \${sample}

    # Handle multiple flowcells per sample
    outSAMattrRGline=""
    readOneFilesIn=""
    readTwoFilesIn=""
    for flowcell in `ls -1 \${sample}/split_\${sample}_*_R1.fastq.gz \
        | xargs -L 1 basename \
        | sed s/split_// \
        | sed s/\${sample}// \
        | cut -d "_" -f 2 \
        | uniq`; do

        # Handle multiple lanes per flowcell
        for lane in `ls -1 \${sample}/split_\${sample}_\${flowcell}_*_R1.fastq.gz \
            | xargs -L 1 basename \
            | sed s/split_// \
            | sed s/\${sample}// \
            | cut -d "_" -f 3 \
            | uniq`; do

            # Construct the read file paths
            if [ -n "\${readOneFilesIn}" ]; then
              readOneFilesIn="\${readOneFilesIn},"
            fi
            readOneFilesIn="\${readOneFilesIn}\${sample}/split_\${sample}_\${flowcell}_\${lane}_R1.fastq.gz"
            if [ -n "\${readTwoFilesIn}" ]; then
              readTwoFilesIn="\${readTwoFilesIn},"
            fi
            readTwoFilesIn="\${readTwoFilesIn}\${sample}/split_\${sample}_\${flowcell}_\${lane}_R2.fastq.gz"

            # Find the multiplex barcode
            barcode=`grep "multiplex_barcode" \${sample}/split_\${sample}_\${flowcell}_\${lane}_R1.report.json \
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
    echo "--genomeDir ${params.genome_dir}/STAR_Indices/Homo_sapiens_and_SARS-CoV-2_ensembl_RL-151" >> \${sample}.log
    echo "--outSAMattrRGline \${outSAMattrRGline}" >> \${sample}.log
    echo "--runThreadN ${params.numberOfThreads}" >> \${sample}.log
    echo "--outFileNamePrefix \${sample}/\${sample}_" >> \${sample}.log
    echo "--readFilesIn \${readFilesIn}" >> \${sample}.log
    STAR \
        --twopassMode Basic \
        --limitBAMsortRAM 65000000000 \
        --genomeDir ${params.genome_dir}/STAR_Indices/Homo_sapiens_and_SARS-CoV-2_ensembl_RL-151 \
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
        --outFileNamePrefix \${sample}/\${sample}_ \
        --readFilesIn \${readFilesIn}
    """
}

params.aligned_multiqc_dir = params.aligned_reads_dir + "/alignment_multiQC_report"

process compileAlignedReadsQC {

    label "COVIRT_fastq_to_alignment"

    publishDir params.aligned_multiqc_dir, mode: "copy"

    input:
    file "*" from star_log_files_ch.collect()

    output:
    file "*" into aligned_multiqc_ch

    """
    multiqc -n alignment_multiqc -f ${params.aligned_reads_dir}
    """
}


process generateStarCountsTable {

    label "COVIRT_fastq_to_alignment"

    publishDir params.aligned_reads_dir, mode: "copy"

    input:
    file "*" from reads_counts_files_ch.collect()

    output:
    file "*" into counts_table_file_ch

    """
    #!/usr/local/bin/Rscript

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
params.AvailableMemoryPerThread = "2600M"
params.NumberOfThreads = 10

process sortAndIndexAlignedReads {

    label "COVIRT_fastq_to_alignment"

    publishDir params.aligned_reads_dir, mode: "copy"

    input:
    file aligned_reads_file from aligned_reads_files_ch

    output:
    file "*/*" into sorted_reads_files_ch

    """
    sample=`echo ${aligned_reads_file} | sed s/_Aligned.sortedByCoord.out.bam//`
    echo Sample \$sample

    mkdir \${sample}

    samtools sort \
        -m ${params.AvailableMemoryPerThread} --threads ${params.NumberOfThreads} \
        -o \${sample}/\${sample}_Aligned.sortedByCoord_sorted.out.bam ${aligned_reads_file}
    samtools index \
        -@ ${params.NumberOfThreads} \${sample}/\${sample}_Aligned.sortedByCoord_sorted.out.bam
    """
}


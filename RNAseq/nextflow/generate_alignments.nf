#!/bin/bash

raw_reads_files_ch = Channel.fromPath(params.raw_reads_dir + "/*.fastq.gz")

params.raw_fastqc_dir = params.COVIRT_Code + "/00-RawData/FastQC_Reports"

process createRawReadsQC {

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=6:model=has -l walltime=5:00:00 -N COVIRT_raw_fastqc -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 6
    //memory = 8.GB

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

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=24:model=has -l walltime=5:00:00 -N COVIRT_raw_multiqc -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 24
    //memory = 16.GB

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

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=10:model=has -l walltime=50:00:00 -N COVIRT_trimmomatic -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 10
    //memory = 40.GB

    label "COVIRT_fastq_to_alignment"
  
    publishDir params.trimmed_reads_dir, mode: "copy", pattern: "*_trimmed.fq.gz"
    publishDir params.trimmed_logs_dir, mode: "copy", pattern: "*_trimming*"

    input:
    set sample, file(raw_reads_file_pair) from raw_reads_file_pairs_ch

    output:
    file "*_P_trimmed.fq.gz" into trimmed_P_reads_files_ch
    file "*_U_trimmed.fq.gz" into trimmed_U_reads_files_ch
    file "*_R1_P_trimmed.fq.gz" into trimmed_reads_one_files_ch
    file "*_R2_P_trimmed.fq.gz" into trimmed_reads_two_files_ch
    file "*_trimming*" into trimming_log_files_ch


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
        ILLUMINACLIP:${params.trimmomatic_adapters}/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
        LEADING:20 \
        TRAILING:20 \
        MINLEN:15
    """
}

params.trimmed_fastqc_dir = params.trimmed_data_dir + "/FastQC_Reports"

params.trimmed_P_fastqc_dir = params.trimmed_fastqc_dir + "/paired_reads"


process createTrimmedPReadsQC {

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=6:model=has -l walltime=5:00:00 -N COVIRT_trimmed_fastqc_paired -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 6
    //memory = 8.GB

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

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=6:model=has -l walltime=5:00:00 -N COVIRT_trimmed_fastqc_unpaired -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 6
    //memory = 8.GB

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

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=24:model=has -l walltime=5:00:00 -N COVIRT_trimmed_multiqc_paired -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 24
    //memory = 16.GB

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

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=24:model=has -l walltime=5:00:00 -N COVIRT_trimmed_multiqc_unpaired -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 24
    //memory = 16.GB

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

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=6:model=has -l walltime=24:00:00 -N COVIRT_fq_split_unfiltered -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 6
    //memory = 28.GB

    // label "COVIRT_fq_splitter"

    publishDir params.split_reads_dir, mode: "copy"

    input:
    file trimmed_reads_one_file from trimmed_reads_one_files_ch
    file trimmed_reads_two_file from trimmed_reads_two_files_ch

    output:
    env sample into split_reads_sample_ch1, split_reads_sample_ch2
    file "*/split_*_R1.fq.gz" into split_reads_one_files_ch
    file "*/split_*_R2.fq.gz" into split_reads_two_files_ch
    file "*/split_*.json" into split_reads_json_files_ch

  """
  sample=`echo ${trimmed_reads_one_file} | sed s/_R1_P_trimmed.fq.gz//`
  echo Sample \$sample

  mkdir \${sample}
  # Activate gdc-fastq-splitter
  module load /nobackupp16/swbuild/hsp/COVID19/anaconda3.modulefile
  cd /nobackupp16/swbuild/hsp/COVID19/gdc-fastq-splitter
  source ./venv/bin/activate
  cd -

  echo ""
  echo "gdc-fastq-splitter version:"
  gdc-fastq-splitter --version
  echo ""

  # source /Users/asaravia/Desktop/Logyx_NASA/GeneLab/COVID19/NASA_HEC_analyses/Pipelines/GitHub_COV-IRT_Pipeline/Conda_Environments/gdc-fastq-splitter/venv/bin/activate
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

  # Deactivate gdc-fastq-splitter
  deactivate
  """
}

params.filtered_reads_dir = params.trimmed_data_dir + "/Fastq_RG_NOrRNA"

params.htstream_logs_dir = params.trimmed_data_dir + "/HTStream_logs"

process filterSplitReads {

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=24:model=has -l walltime=24:00:00 -N COVIRT_HTStream_rRNA_fqsplit -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 24
    //memory = 60.GB

    label "COVIRT_HTStream"

    publishDir params.filtered_reads_dir, mode: "copy", pattern: "*/*.fastq.gz"
    publishDir params.htstream_logs_dir, mode: "copy", pattern: "*_htsStats.log"

    input:
    env sample from split_reads_sample_ch1
    file split_reads_one_file from split_reads_one_files_ch
    file split_reads_two_file from split_reads_two_files_ch

    output:
    env sample into filtered_reads_sample_ch
    file "*/split_*_R1.fastq.gz" into filtered_reads_one_files_ch
    file "*/split_*_R2.fastq.gz" into filtered_reads_two_files_ch
    file "split_*.log" into filtered_reads_log_files_ch

    """
    mkdir \${sample}

    # Handle multiple flowcells per sample
    for flowcell in `ls -1 split_\${sample}_*_R1.fq.gz \
        | xargs -L 1 basename \
        | sed s/split_// \
        | sed s/\${sample}// \
        | cut -d "_" -f 2 \
        | uniq`; do

        # Handle multiple lanes per flowcell
        for lane in `ls -1 split_\${sample}_\${flowcell}_*_R1.fq.gz \
            | xargs -L 1 basename \
            | sed s/split_// \
            | sed s/\${sample}// \
            | cut -d "_" -f 3 \
            | uniq`; do
            reads_one_file="split_\${sample}_\${flowcell}_\${lane}_R1.fq.gz"
            reads_two_file="split_\${sample}_\${flowcell}_\${lane}_R2.fq.gz"
            base_name="split_\${sample}_\${flowcell}_\${lane}"
            # TODO: Need to get the reference file to the worker
            hts_SeqScreener \
                -L \${base_name}_htsStats.log \
                -1 \${reads_one_file} \
                -2 \${reads_two_file} \
                -s ${params.rRNA_ref} \
                -x 0.20 \
                -f \${sample}/\${base_name}
        done
    done
    """
}


// Generate alignment data without rRNA - removed in silico

params.NOrRNA_aligned_reads_dir = params.COVIRT_Code + "/02-AlignedData_NOrRNA"

process NOrRNAalignReadsToReferenceGenome {

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=24:model=has -l walltime=72:00:00 -N COVIRT_STAR_alignment -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 24
    //memory = 80.GB

    label "COVIRT_fastq_to_align_STAR_2_7_1a"

    publishDir params.NOrRNA_aligned_reads_dir, mode: "copy"

    input:
    env sample from filtered_reads_sample_ch
    //file split_reads_json_file from split_reads_json_files_ch

    output:
    file "*/*ReadsPerGene.out.tab" into NOrRNA_reads_counts_files_ch
    file "*/*_Aligned.sortedByCoord.out.bam" into NOrRNA_aligned_reads_files_ch
    file "*/*_Aligned.toTranscriptome.out.bam" into NOrRNA_aligned_transcriptome_reads_files_ch
    file "*/*_Log.final.out" into NOrRNA_star_log_files_ch
    file "*/*_Chimeric.out.junction" into NOrRNA_chimeric_junct_ch
    file "*/*_Chimeric.out.sam" into NOrRNA_chimeric_sam_ch
    file "*/*_SJ.out.tab" into NOrRNA_SJ_table_ch

    """
    mkdir \${sample}

    # Copy in reads files (can't use links anyway)
    cp ${params.filtered_reads_dir}/\${sample}/split_\${sample}_*.fastq.gz \${sample}
    cp ${params.split_reads_dir}/\${sample}/split_\${sample}_*.report.json \${sample}

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
    echo "--genomeDir ${params.STAR_index}" >> \${sample}.log
    echo "--outSAMattrRGline \${outSAMattrRGline}" >> \${sample}.log
    echo "--runThreadN ${params.numberOfThreads}" >> \${sample}.log
    echo "--outFileNamePrefix \${sample}/\${sample}_" >> \${sample}.log
    echo "--readFilesIn \${readFilesIn}" >> \${sample}.log
    STAR \
        --twopassMode Basic \
        --limitBAMsortRAM 65000000000 \
        --genomeDir ${params.STAR_index} \
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

params.NOrRNA_aligned_multiqc_dir = params.NOrRNA_aligned_reads_dir + "/alignment_multiQC_report"

process NOrRNAcompileAlignedReadsQC {

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=24:model=has -l walltime=10:00:00 -N COVIRT_align_multiqc -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 24
    //memory = 18.GB

    label "COVIRT_fastq_to_alignment"

    publishDir params.NOrRNA_aligned_multiqc_dir, mode: "copy"

    input:
    file "*" from NOrRNA_star_log_files_ch.collect()

    output:
    file "*" into NOrRNA_aligned_multiqc_ch

    """
    multiqc -n alignment_multiqc -f ${params.NOrRNA_aligned_reads_dir}
    """
}


process NOrRNAgenerateStarCountsTable {

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=24:model=has -l walltime=18:00:00 -N COVIRT_STAR_rawCounts_table -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 24
    //memory = 20.GB

    label "COVIRT_fastq_to_alignment"

    publishDir params.NOrRNA_aligned_reads_dir, mode: "copy"

    input:
    file "*" from NOrRNA_reads_counts_files_ch.collect()

    output:
    file "*" into NOrRNA_counts_table_file_ch

    """
    #!/nobackupp16/swbuild/hsp/COVID19/R-3.6.0/bin/Rscript --no-save

    R.Version()

    print("Make STAR counts table")
    print("")

    ### Import data
    ##ff <- list.files(
    ##    file.path("${params.NOrRNA_aligned_reads_dir}"),
    ##    pattern="*/*ReadsPerGene.out.tab", full.names=TRUE
    ##)

    ### Import data
    ff <- list.files(
        pattern="ReadsPerGene.out.tab", full.names=TRUE
    )

    print("check ff")
    ff

    ## Remove the first 4 lines
    counts.files <- lapply(ff, read.table, skip = 4)

    print("check counts.files")
    head(counts.files)

    ## Get counts aligned to the second, reverse, strand
    counts <- as.data.frame(sapply(counts.files, function(x) x[ , 4 ]))

    print("check counts")
    head(counts)

    ## Get sample names
    samples <- sub(
        "_ReadsPerGene.out.tab", "", basename(ff)
    )

    print("check samples")
    samples

    ## Add column and row names
    colnames(counts) <- samples
    row.names(counts) <- counts.files[[1]]\$V1

    print("check counts with col and row labels")
    head(counts)

    ### Export unnormalized counts table
    write.csv(counts, file='STAR_Unnormalized_Counts.csv')

    print("Session Info below: ")
    print("")
    sessionInfo()
    quit()
    """
}

// TODO: Automate setting of these values
params.AvailableMemoryPerThread = "2600M"
params.NumberOfThreads = 10

process NOrRNAsortAndIndexAlignedReads {

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=12:model=has -l walltime=24:00:00 -N COVIRT_samtools_sort -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 12
    //memory = 40.GB

    label "COVIRT_fastq_to_alignment"

    publishDir params.NOrRNA_aligned_reads_dir, mode: "copy"

    input:
    file aligned_reads_file from NOrRNA_aligned_reads_files_ch

    output:
    file "*/*" into NOrRNA_sorted_reads_files_ch

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


// Generate alignment data with rRNA - no in silico rRNA removal

params.aligned_reads_dir = params.COVIRT_Code + "/02-AlignedData"

process alignReadsToReferenceGenome {

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=24:model=has -l walltime=72:00:00 -N COVIRT_STAR_alignment -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 24
    //memory = 80.GB

    label "COVIRT_fastq_to_align_STAR_2_7_1a"

    publishDir params.aligned_reads_dir, mode: "copy"

    input:
    env sample from split_reads_sample_ch2
    //file split_reads_json_file from split_reads_json_files_ch

    output:
    file "*/*ReadsPerGene.out.tab" into reads_counts_files_ch
    file "*/*_Aligned.sortedByCoord.out.bam" into aligned_reads_files_ch
    file "*/*_Aligned.toTranscriptome.out.bam" into aligned_transcriptome_reads_files_ch
    file "*/*_Log.final.out" into star_log_files_ch
    file "*/*_Chimeric.out.junction" into chimeric_junct_ch
    file "*/*_Chimeric.out.sam" into chimeric_sam_ch
    file "*/*_SJ.out.tab" into SJ_table_ch

    """
    mkdir \${sample}

    # Copy in reads files (can't use links anyway)
    cp ${params.split_reads_dir}/\${sample}/split_\${sample}_*.fq.gz \${sample}
    cp ${params.split_reads_dir}/\${sample}/split_\${sample}_*.report.json \${sample}

    # Handle multiple flowcells per sample
    outSAMattrRGline=""
    readOneFilesIn=""
    readTwoFilesIn=""
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

            # Construct the read file paths
            if [ -n "\${readOneFilesIn}" ]; then
              readOneFilesIn="\${readOneFilesIn},"
            fi
            readOneFilesIn="\${readOneFilesIn}\${sample}/split_\${sample}_\${flowcell}_\${lane}_R1.fq.gz"
            if [ -n "\${readTwoFilesIn}" ]; then
              readTwoFilesIn="\${readTwoFilesIn},"
            fi
            readTwoFilesIn="\${readTwoFilesIn}\${sample}/split_\${sample}_\${flowcell}_\${lane}_R2.fq.gz"

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
    echo "--genomeDir ${params.STAR_index}" >> \${sample}.log
    echo "--outSAMattrRGline \${outSAMattrRGline}" >> \${sample}.log
    echo "--runThreadN ${params.numberOfThreads}" >> \${sample}.log
    echo "--outFileNamePrefix \${sample}/\${sample}_" >> \${sample}.log
    echo "--readFilesIn \${readFilesIn}" >> \${sample}.log
    STAR \
        --twopassMode Basic \
        --limitBAMsortRAM 65000000000 \
        --genomeDir ${params.STAR_index} \
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

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=24:model=has -l walltime=10:00:00 -N COVIRT_align_multiqc -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 24
    //memory = 18.GB

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

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=24:model=has -l walltime=18:00:00 -N COVIRT_STAR_rawCounts_table -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 24
    //memory = 20.GB

    label "COVIRT_fastq_to_alignment"

    publishDir params.aligned_reads_dir, mode: "copy"

    input:
    file "*" from reads_counts_files_ch.collect()

    output:
    file "*" into counts_table_file_ch

    """
    #!/nobackupp16/swbuild/hsp/COVID19/R-3.6.0/bin/Rscript --no-save

    R.Version()

    print("Make STAR counts table")
    print("")

    ### Import data
    ##ff <- list.files(
    ##    file.path("${params.aligned_reads_dir}"),
    ##    pattern="*/*ReadsPerGene.out.tab", full.names=TRUE
    ##)

    ### Import data
    ff <- list.files(
        pattern="ReadsPerGene.out.tab", full.names=TRUE
    )

    print("check ff")
    ff

    ## Remove the first 4 lines
    counts.files <- lapply(ff, read.table, skip = 4)

    print("check counts.files")
    head(counts.files)

    ## Get counts aligned to the second, reverse, strand
    counts <- as.data.frame(sapply(counts.files, function(x) x[ , 4 ]))

    print("check counts")
    head(counts)

    ## Get sample names
    samples <- sub(
        "_ReadsPerGene.out.tab", "", basename(ff)
    )

    print("check samples")
    samples

    ## Add column and row names
    colnames(counts) <- samples
    row.names(counts) <- counts.files[[1]]\$V1

    print("check counts with col and row labels")
    head(counts)

    ### Export unnormalized counts table
    write.csv(counts, file='STAR_Unnormalized_Counts.csv')

    print("Session Info below: ")
    print("")
    sessionInfo()
    quit()
    """
}

// TODO: Automate setting of these values
params.AvailableMemoryPerThread = "2600M"
params.NumberOfThreads = 10

process sortAndIndexAlignedReads {

    // Set PBS options
    executor = 'pbspro'
    clusterOptions = '-W group_list=e2255 -S /bin/bash -lselect=1:ncpus=12:model=has -l walltime=24:00:00 -N COVIRT_samtools_sort -q covid19 -M amanda.m.saravia-butler@nasa.gov -m abe -j oe'
    //cpus = 12
    //memory = 40.GB

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


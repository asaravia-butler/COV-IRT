#!/usr/bin/env nextflow

// TODO: Update based on NASA storage conventions
params.aligned_bams_dir = "/data/home/snagar9/data/covirt-nextflow/data/Fastq_Input_Files_for_Testing/aligned_reads/"

// TODO: Integrate into worflow so that the Channel doesn't need to be populated from a path
aligned_reads_files_ch = Channel.fromPath(params.aligned_bams_dir + "*_Aligned.sortedByCoord_sorted.out.bam")

// TODO: Change paths as needed
params.variant_calling_op_dir = "/data/home/snagar9/data/covirt-nextflow/data/variant_calling"

// TODO: Automate setting of this value
params.numberOfThreads = 16

process markDuplicates {
    // TODO: Uncomment
    // label: "covirt_gatk"
    
    publishDir params.variant_calling_op_dir, mode: "copy"

    input:
        file aligned_reads_file from aligned_reads_files_ch
    
    output:
        file "*marked_duplicates.bam" into dupe_marked_bams_ch
        file "*marked_dup_metrics.txt" into dupes_marked_metrics_ch
        file "*marked_duplicates.bam.bai" into dupe_marked_bams_index_ch

    """
    sample=`echo ${aligned_reads_file} | sed 's/_Aligned.sortedByCoord_sorted.out.bam//'`

    echo Sample \$sample

    # Marking duplicate reads
    ~/setups/gatk-4.0.10.1/gatk MarkDuplicates -I ${aligned_reads} \
            -O ./\${sample}_marked_duplicates.bam \
            -M ./\${sample}_marked_dup_metrics.txt

    # Indexing duplicate reads
    samtools index -@ ${params.numberOfThreads} ./\${sample}_marked_duplicates.bam
    """
}

// TODO: Write possible control flow to choose either human-filtered reads or unfiltered data
params.ref_genome = "~/data/covirt-nextflow/data/Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa"

process splitReads {
    // TODO: Uncomment
    // label: "covirt_gatk"

    publishDir params.variant_calling_op_dir, mode: "copy"

    input:
        file dupe_marked_bam from dupe_marked_bams_ch
    
    output:
        file "*Split.bam" into cigar_split_bams_ch
    
    """
    sample=`echo ${aligned_reads_file} | sed 's/_marked_duplicates.bam//'`

    echo Sample \$sample

    ~/setups/gatk-4.0.10.1/gatk --java-options "-Xmx100G" SplitNCigarReads -R ${params.ref_genome} \
        -I ${dupe_marked_bam} \
        -O ./\${sample}_Split.bam
    """

}

// TODO: Change path
params.indel_ref_loc = "~/data/covirt-nextflow/data/Homo_sapiens_assembly38.ens100.known_indels.vcf.gz"
params.dbsnp_loc = "~/data/covirt-nextflow/data/dbSNP_v153_ens.vcf.gz"

// TODO: Possible control flow for human-filtered or unfiltered analysis
params.ref_genome_dict = "~/data/covirt-nextflow/data/Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.dict"

process generateRecalTable {
    // TODO: Uncomment
    // label: "covirt_gatk"

    publishDir params.variant_calling_op_dir, mode: "copy"

    input:
        file cigar_split_bam from cigar_split_bams_ch
    
    output:
        file "*recal_data.table" into recal_tables_ch

    """
    sample=`echo ${aligned_reads_file} | sed 's/_Split.bam//'`

    echo Sample \$sample

    ~/setups/gatk-4.0.10.1/gatk --java-options "-Xmx100G" BaseRecalibrator -I ${cigar_split_bam} \
        -R ${params.ref_genome} \
        --known-sites ${params.indel_ref_loc} \
        --known-sites ${params.dbsnp_loc} \
        -O ./\${sample}_recal_data.table \
        --sequence-dictionary ${params.ref_genome_dict}
    """
}

# Sorry for the long name... \o/
process compareBaseQualityScoreRecalTables {
    // TODO: Uncomment
    // label: "covirt_gatk"

    publishDir params.variant_calling_op_dir, mode: "copy"

    input:
        file recal_table from recal_tables_ch
    
    output:
        file "*AnalyzeCovariates.pdf" into comparison_reports_ch
    
    """
    sample=`echo ${aligned_reads_file} | sed 's/_Split.bam//'`

    echo Sample \$sample

    ~/setups/gatk-4.0.10.1/gatk --java-options "-Xmx100G" AnalyzeCovariates -bqsr ${recal_table} \
          -plots ./\${sample}_AnalyzeCovariates.pdf
    """
}

process applyBQSR {
    // TODO: Uncomment
    // label: "covirt_gatk"

    publishDir params.variant_calling_op_dir, mode: "copy"

    input:
        file recal_table from recal_tables_ch
        file cigar_split_bam from cigar_split_bams_ch
    
    output:
        file "*BSQR-applied.out.bam" into bqsr_bams_ch
    
    """
    sample=`echo ${aligned_reads_file} | sed 's/_Split.bam//'`

    echo Sample \$sample

    gatk --java-options "-Xmx100G" ApplyBQSR -R ${params.ref_genome} \
        -I ${cigar_split_bam} \
        --bqsr-recal-file ${recal_table} \
        -O ./\${sample}_BSQR-applied.out.bam

}

process callVariants {
    // TODO: Uncomment
    // label: "covirt_gatk"

    publishDir params.variant_calling_op_dir, mode: "copy"

    input:
        // TODO: Check whether we're okay doing only autosomes
        each chr from 1..22
        file bqsr_bam from bqsr_bams_ch
    
    output:
        file "*vcf.gz" into gvcf_ch
    """
    sample=`echo ${aligned_reads_file} | sed 's/_BSQR-applied.out.bam//'`

    echo Sample \$sample

    gatk --java-options "-Xmx100G" HaplotypeCaller -R ${params.ref_genome} \
        -I ${bqsr_bam} \
        -O ./\${sample}_${chr}.vcf.gz \
        -ERC GVCF \
        --dont-use-soft-clipped-bases \
        --standard-min-confidence-threshold-for-calling 20 \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --intervals ${chr}
    """
}

process makeGenomicsDB {
    // TODO: Uncomment
    // label: "covirt_gatk"

    publishDir params.variant_calling_op_dir, mode: "copy"

    input:
        // TODO: Check whether we're okay doing only autosomes
        each chr in 1..22
        file all_gvcfs from gvcf_ch.collect()
    
    output:
        file "*vcf.gz" into gvcf_ch
    
    """
    # Generating sample map of all generated VCF files
    

}

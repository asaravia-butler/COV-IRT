#!/usr/bin/env nextflow

params.aligned_bams_dir = params.COVIRT_Code + "/02-AlignedData"


// TODO: Integrate into worflow so that the Channel doesn't need to be populated from a path
aligned_reads_files_ch = Channel.fromPath(params.aligned_bams_dir + "/*/*_Aligned.sortedByCoord_sorted.out.bam")

params.variant_calling_dir = params.COVIRT_Code + "/GATK_analysis"

params.variant_calling_intermed_dir = params.variant_calling_dir + "/03-GATKIntermediates/"
params.variant_calling_output_dir = params.variant_calling_dir + "/04-GATKOutput/"

// TODO: Automate setting of this value
params.numberOfThreads = 10

process markDuplicates {

    label "COVIRT_GATK"
    
    publishDir params.variant_calling_intermed_dir, mode: "copy"

    input:
    file aligned_reads_file from aligned_reads_files_ch
    
    output:
    file "*/*marked_duplicates.bam" into dupe_marked_bams_ch
    file "*/*marked_dup_metrics.txt" into dupes_marked_metrics_ch
    file "*/*marked_duplicates.bam.bai" into dupe_marked_bams_index_ch

    """
    sample=`echo ${aligned_reads_file} | sed 's/_Aligned.sortedByCoord_sorted.out.bam//'`
    echo Sample \$sample

    mkdir \${sample}

    # Marking duplicate reads
    gatk MarkDuplicates \
        -I ${aligned_reads_file} \
        -O \${sample}/\${sample}_marked_duplicates.bam \
        -M \${sample}/\${sample}_marked_dup_metrics.txt

    # Indexing duplicate reads
    samtools index -@ ${params.numberOfThreads} \${sample}/\${sample}_marked_duplicates.bam
    """
}


// TODO: Write possible control flow to choose either human-filtered reads or unfiltered data
// params.ref_genome = "~/data/covirt-nextflow/data/Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa"
params.ref_genome = params.COVIRT_Data + "/Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa"

//dupe_marked_bams_ch = Channel.fromPath(params.aligned_bams_dir + "*/*_marked_duplicates.bam")

process splitReads {

    label "COVIRT_GATK"

    publishDir params.variant_calling_intermed_dir, mode: "copy"

    input:
    file dupe_marked_bam from dupe_marked_bams_ch
    
    output:
    file "*/*Split.bam" into cigar_split_bams_1_ch
    
    """
    sample=`echo ${dupe_marked_bam} | sed 's/_marked_duplicates.bam//'`
    echo Sample \$sample

    mkdir \${sample}

    gatk --java-options "-Xmx100G" SplitNCigarReads \
        -R ${params.ref_genome} \
        -I ${dupe_marked_bam} \
        -O \${sample}/\${sample}_Split.bam
    """
}

params.indel_ref_loc = params.COVIRT_Data + "/Homo_sapiens_assembly38.ens100.known_indels.vcf.gz"
params.dbsnp_loc = params.COVIRT_Data + "/dbSNP_v153_ens.vcf.gz"

// TODO: Possible control flow for human-filtered or unfiltered analysis
// params.ref_genome_dict = "~/data/covirt-nextflow/data/Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.dict"
params.ref_genome_dict = params.COVIRT_Data + "/Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.dict"

process generateRecalTable {

    label "COVIRT_GATK"

    publishDir params.variant_calling_intermed_dir, mode: "copy", pattern: "*/*recal_data.table"

    input:
    file cigar_split_bam from cigar_split_bams_1_ch
    
    output:
    file "*/*recal_data.table" into recal_tables_1_ch, recal_tables_2_ch
    file cigar_split_bam into cigar_split_bams_2_ch

    """
    sample=`echo ${cigar_split_bam} | sed 's/_Split.bam//'`
    echo Sample \$sample

    mkdir \${sample}

    gatk --java-options "-Xmx100G" BaseRecalibrator \
        -I ${cigar_split_bam} \
        -R ${params.ref_genome} \
        --known-sites ${params.indel_ref_loc} \
        --known-sites ${params.dbsnp_loc} \
        -O \${sample}/\${sample}_recal_data.table \
        --sequence-dictionary ${params.ref_genome_dict}
    """
}

process compareBaseQualityScoreRecalTables {

    label "COVIRT_GATK"

    publishDir params.variant_calling_intermed_dir, mode: "copy"

    input:
    file recal_table from recal_tables_1_ch
    
    output:
    file "*/*AnalyzeCovariates.pdf" into comparison_reports_ch
    
    """
    sample=`echo ${recal_table} | sed 's/_recal_data.table//'`
    echo Sample \$sample

    mkdir \${sample}

    gatk --java-options '-Xmx100G' AnalyzeCovariates \
        -bqsr ${recal_table} \
        -plots \${sample}/\${sample}_AnalyzeCovariates.pdf
    """
}

process applyBQSR {

    label "COVIRT_GATK"

    publishDir params.variant_calling_intermed_dir, mode: "copy"

    input:
    file recal_table from recal_tables_2_ch
    file cigar_split_bam from cigar_split_bams_2_ch
    
    output:
    file "*/*BSQR-applied.out.bam" into bqsr_bams_ch
    
    """
    sample=`echo ${recal_table} | sed 's/_recal_data.table//'`
    echo Sample \$sample

    mkdir \${sample}

    gatk --java-options "-Xmx100G" ApplyBQSR \
        -R ${params.ref_genome} \
        -I ${cigar_split_bam} \
        --bqsr-recal-file ${recal_table} \
        -O \${sample}/\${sample}_BSQR-applied.out.bam
    """

}

chrs_1_ch = Channel.from(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                         21, 22, "X", "Y", "MT")

//chrs_1_ch.view()

process callVariants {

    label "COVIRT_GATK"

    publishDir params.variant_calling_output_dir, mode: "copy"

    input:
    // TODO: Check whether we're okay doing only autosomes
    each chr from chrs_1_ch
    file bqsr_bam from bqsr_bams_ch
    
    output:
    file "*/*vcf.gz" into gvcf_ch
    file "*/*vcf.gz.tbi" into gvcf_tbi_ch

    """
    sample=`echo ${bqsr_bam} | sed 's/_BSQR-applied.out.bam//'`
    echo Sample \$sample

    mkdir \${sample}

    gatk --java-options "-Xmx100G" HaplotypeCaller \
        -R ${params.ref_genome} \
        -I ${bqsr_bam} \
        -O \${sample}/\${sample}_chr${chr}.vcf.gz \
        -ERC GVCF \
        --dont-use-soft-clipped-bases \
        --standard-min-confidence-threshold-for-calling 20 \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --intervals ${chr}
    """
}

params.genomics_db = params.variant_calling_output_dir + "/GVCF_databases/"

chrs_2_ch = Channel.from(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                         21, 22, "X", "Y", "MT")

//chrs_2_ch.view()

process makeGenomicsDB {

    label "COVIRT_GATK"

    publishDir params.genomics_db, mode: "copy"

    input:
    // TODO: Check whether we're okay doing only autosomes
    each chr from chrs_2_ch
    file all_gvcfs from gvcf_ch.collect()
    file all_gvcf_tbis from gvcf_tbi_ch.collect()
    
    output:
    file "*database" into genomics_db_ch
    
    """
    # Generate sample map of all generated VCF files
    echo ${all_gvcfs} \
        | tr ' ' '\n' \
        | grep _chr${chr}.vcf.gz > chr${chr}_column2.txt

    echo ${all_gvcfs} \
        | tr ' ' '\n' \
        | grep _chr${chr}.vcf.gz \
        | awk -F_chr '{sub("_chr"\$NF,""); print}' > chr${chr}_column1.txt

    paste -d"\t" chr${chr}_column1.txt chr${chr}_column2.txt > sample_map_${chr}.txt
    
    # TODO: Remove stacketrace on user exception
    gatk --java-options "-Xmx40G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport \
        --sample-name-map sample_map_${chr}.txt \
        --genomicsdb-workspace-path chr${chr}_database \
        --intervals ${chr}
    """
}

process jointGenotyping {

    label "COVIRT_GATK"

    publishDir params.variant_calling_output_dir + "/all_samples", mode: "copy"

    input:
    file genomics_db from genomics_db_ch

    output:
    file "*/*Geno_out.vcf.gz" into joint_called_vcf_ch
    file "*/*Geno_out.vcf.gz.tbi" into joint_called_tbi_ch

    """
    # Getting chromosome number
    chr_num=`echo ${genomics_db} | sed 's/_database//; s/chr//'`

    mkdir chr\${chr_num}

    gatk --java-options "-Xmx40G" GenotypeGVCFs \
       -R ${params.ref_genome} \
       -V gendb://${genomics_db} \
       -G StandardAnnotation \
       -G AS_StandardAnnotation \
       -O chr\${chr_num}/chr\${chr_num}_Geno_out.vcf.gz
    """
}

process variantAnnotationFilter {

    label "COVIRT_GATK"

    publishDir params.variant_calling_output_dir + "/all_samples", mode: "copy"

    input:
    file joint_called_vcf from joint_called_vcf_ch
    file joint_called_tbi from joint_called_tbi_ch

    output:
    file "*/*VarFilt_output.vcf.gz" into annot_filtered_vcf_ch
    file "*/*VarFilt_output.vcf.gz.tbi" into annot_filtered_tbi_ch

    """
    chr_num=`echo ${joint_called_vcf} | sed 's/_Geno_out\\.vcf\\.gz//; s/chr//'`

    mkdir chr\${chr_num}

    gatk VariantFiltration -R ${params.ref_genome} \
        -V ${joint_called_vcf} \
        -O chr\${chr_num}/chr\${chr_num}_VarFilt_output.vcf.gz \
        --window 35 \
        --cluster 3 \
        --filter-name "FS" \
        --filter "FS > 30.0" \
        --filter-name "QD" \
        --filter "QD < 2.0"
    """
}

process combineVCFs {

    label "COVIRT_GATK"

    publishDir params.variant_calling_output_dir + "/all_samples", mode: "copy"

    input:
    file all_filt_vcfs from annot_filtered_vcf_ch.collect()
    file all_filt_tbis from annot_filtered_tbi_ch.collect()

    output:
    file "merged_chr.vcf.gz" into merged_vcf_ch

    """
    # Create the input file list
    for gvcf in ${all_filt_vcfs}; do
        echo \${gvcf} >> filt_vcfs.list
    done
    gatk MergeVcfs \
        --INPUT filt_vcfs.list \
        --SEQUENCE_DICTIONARY ${params.ref_genome_dict} \
        --OUTPUT merged_chr.vcf.gz
    """
}


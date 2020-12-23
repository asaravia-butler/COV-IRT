# GATK Variant Calls Processing Pipeline

> **This page holds an overview and instructions for how COV-IRT generates variant call data from raw Illumina RNA-sequencing data of COVID-19 samples. This pipeline uses the indexed \*Aligned.sortedByCoord_sorted.out.bam files generated from [step 6 of the 'Raw to Aligned Data Pipeline'](Raw_to_Aligned_Data_Pipeline.md#6-sort-and-index-genome-aligned-data) and the samples.txt file generated from [step 0 of the 'Raw to Aligned Data Pipeline'](Raw_to_Aligned_Data_Pipeline.md#0-create-sample-list). Exact processing commands used for specific cohorts of samples are available in the [Exact_scripts_used](Exact_scripts_used) sub-directory.**  

---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Retrieve and Index Genome Files and Create Reference Dictionary**](#1-retrieve-and-index-genome-files-and-create-reference-dictionary)
    - [**1a. Get Genome Files**](#1a-get-genome-files)
    - [**1b. Index Genome Files**](#1b-index-genome-files)
    - [**1c. Create Reference Dictionaries**](#1c-create-reference-dictionaries)
  - [**2. Retrieve Variant Call Reference Files and Convert to Match Genome IDs**](#2-retrieve-variant-call-reference-files-and-convert-to-match-genome-ids)
    - [**2a. Get Variant Call Reference Files**](#2a-get-variant-call-reference-files)
    - [**2b. Convert Variant Call Reference Files to Match Genome IDs then Index**](#2b-convert-variant-call-reference-files-to-match-genome-ids-then-index)
  - [**3. Define Paths to Directories Containing Input, Reference, Intermediate, and Output Files**](#3-define-paths-to-directories-containing-input-reference-intermediate-and-output-files)
  - [**4. Mark and Index Duplicate Reads**](#4-mark-and-index-duplicate-reads)
    - [**4a. Mark Duplicate Reads**](#4a-mark-duplicate-reads)
    - [**4b. Index Duplicate Reads**](#4b-index-duplicate-reads)
  - [**5. Split Reads with N in Cigar**](#5-split-reads-with-n-in-cigar)
  - [**6. Generate Recalibration Table for Base Quality Score Recalibration**](#6-generate-recalibration-table-for-base-quality-score-recalibration)
  - [**7. Evaluate and Compare Base Quality Score Recalibration Tables**](#7-evaluate-and-compare-base-quality-score-recalibration-tables)
  - [**8. Apply Base Quality Score Recalibration**](#8-apply-base-quality-score-recalibration)
  - [**9. Call Germline SNPs and Indels via Local Re-assembly of Haplotypes**](#9-call-germline-snps-and-indels-via-local-re-assembly-of-haplotypes)
  - [**10. Import VCFs to Genomics Database**](#10-import-vcfs-to-genomics-database)
  - [**11. Perform Joint Genotyping**](#11-perform-joint-genotyping)
  - [**12. Filter Variant Calls Based on Annotations**](#12-filter-variant-calls-based-on-annotations)
  - [**13. Combine Variant Files From Each Chromosome**](#13-combine-variant-files-from-each-chromosome)

---

<br>

# Software used  

|Program|Version*|Relevant Links|
|:------|:------:|:-------------|
|Samtools|`samtools --version`|[http://www.htslib.org/](http://www.htslib.org/)|
|gatk|`gatk --version`|[https://gatk.broadinstitute.org/](https://gatk.broadinstitute.org/)|
|fgbio|`fgbio --version`|[http://fulcrumgenomics.github.io/fgbio/](http://fulcrumgenomics.github.io/fgbio/)|

>**\*** Exact versions used to process specific cohorts are available in the [Exact_scripts_used](Exact_scripts_used) sub-directory. 

---

<br>

# General processing overview with example commands  

> Exact processing commands used for specific cohorts are provided in the [Exact_scripts_used](Exact_scripts_used) sub-directory.  

---

## 1. Retrieve and Index Genome Files and Create Reference Dictionary

### 1a. Get Genome Files 

> Genome files used in this pipeline are the same as those used in [step 4 of the 'Raw to Aligned Data Pipeline'](Raw_to_Aligned_Data_Pipeline.md#4-retrieve-genomeannotation-files-and-build-star-reference).

Get human fasta file from [Ensembl](https://www.ensembl.org/) - used for processing human-filtered data and needed to process unfiltered data

```
wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

Get SARS-CoV-2 fasta file from [Ensembl](https://www.ensembl.org/) - needed to process unfiltered data

```
wget ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz  
```

Concatenate human and SARS-CoV-2 fasta files - concatenated files used for processing unfiltered data

```
zcat Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa
```

<br>

### 1b. Index Genome Files  

```
samtools faidx /path/to/genome/file
```

**Input Data:**

For processing human-filtered data
- Homo_sapiens.GRCh38.dna.primary_assembly.fa (genome sequence)

For processing unfiltered data
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa (genome sequence)

**Output Data:**

For processing human-filtered data
- Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai (genome index)

For processing unfiltered data
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa.fai (genome index)

<br>

### 1c. Create Reference Dictionaries  

```
gatk --java-options "-Xmx100G" CreateSequenceDictionary \
  -R /path/to/genome/file \
  -O /path/to/genome/index/file
```

**Input Data:**

For processing human-filtered data
- Homo_sapiens.GRCh38.dna.primary_assembly.fa (genome sequence)

For processing unfiltered data
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa (genome sequence)

**Output Data:**

For processing human-filtered data
- Homo_sapiens.GRCh38.dna.primary_assembly.dict (genome dictionary)

For processing unfiltered data
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.dict (genome dictionary)

---

<br>

## 2. Retrieve Variant Call Reference Files and Convert to Match Genome IDs

### 2a. Get Variant Call Reference Files  

Get human known indels file from the GATK resource bundle 

```
wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
```

Get human SNP database from [NCBI](https://www.ncbi.nlm.nih.gov/snp/)

```
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz 
```

<br>

### 2b. Convert Variant Call Reference Files to Match Genome IDs then Index 

Convert human known indels file to match Ensembl genome IDs

```
java -jar fgbio-1.2.0.jar UpdateVcfContigNames -i Homo_sapiens_assembly38.known_indels.vcf.gz \
  --skip-missing true \
  -d ENS100_MN908947.3.dna.primary_assembly.dict \
  -o Homo_sapiens_assembly38.ens100.known_indels.vcf.gz
```

```
tabix -p vcf Homo_sapiens_assembly38.ens100.known_indels.vcf.gz
```

Convert human SNP database to match Ensembl genome IDs

```
java -jar fgbio-1.2.0.jar UpdateVcfContigNames -i 00-All.vcf.gz \
  --skip-missing true \
  -d ENS100_MN908947.3.dna.primary_assembly.dict \
  -o dbSNP_v153_ens.vcf.gz
```

```
tabix -p vcf dbSNP_v153_ens.vcf.gz
```

**Input Data:**

For human known indels 
- Homo_sapiens_assembly38.known_indels.vcf.gz (human indels reference)

For human SNP database
- 00-All.vcf.gz (human SNP database)

**Output Data:**

For human known indels 
- Homo_sapiens_assembly38.ens100.known_indels.vcf.gz (human indels reference with ensembl IDs)
- Homo_sapiens_assembly38.ens100.known_indels.vcf.gz.tbi (human indels index with ensembl IDs)

For human SNP database
- dbSNP_v153_ens.vcf.gz (human SNP database with ensembl IDs)
- dbSNP_v153_ens.vcf.gz.tbi (human SNP database index with ensembl IDs)

---

<br>

## 3. Define Paths to Directories Containing Input, Reference, Intermediate, and Output Files

```
in_dir=/path/to/directory/containing/*Aligned.sortedByCoord_sorted.out.bam/files
inter_dir=/path/to/GATK/intermediate/output/files
out_dir=/path/to/GATK/vcf/output/files
genome_ref=/path/to/directory/containing/indexed_genome/and/genome_dict/files
vcf_ref=/path/to/directory/containing/indexed_variant_call_reference/files
```

---

<br>

## 4. Mark and Index Duplicate Reads

### 4a. Mark Duplicate Reads 

```
gatk MarkDuplicates -I $in_dir/${sample}/${sample}_Aligned.sortedByCoord_sorted.out.bam \
  -O $inter_dir/${sample}/${sample}_marked_duplicates.bam \
  -M $inter_dir/${sample}/${sample}_marked_dup_metrics.txt
```

**Input Data:**
- *Aligned.sortedByCoord_sorted.out.bam (samtools sorted mapping to genome file, output from [step 6a of the 'Raw to Aligned Data Pipeline'](Raw_to_Aligned_Data_Pipeline.md#6a-sort-genome-aligned-data))

**Output Data:**
- *marked_duplicates.bam (samtools sorted mapping to genome with duplicates marked file)
- *marked_dup_metrics.txt (file containing duplication metrics)

<br>

### 4b. Index Duplicate Reads  

```
samtools index -@ NumberOfThreads $inter_dir/${sample}/${sample}_marked_duplicates.bam
```

**Input Data:** 
- *marked_duplicates.bam (samtools sorted mapping to genome with duplicates marked file, output from step 4a)

**Output Data:**
- *marked_duplicates.bam.bai (index of samtools sorted mapping to genome with duplicates marked file)

---

<br>

## 5. Split Reads with N in Cigar

```
gatk --java-options "-Xmx100G" SplitNCigarReads -R $genome_ref/*.fa \
  -I $inter_dir/${sample}/${sample}_marked_duplicates.bam \
  -O $inter_dir/${sample}/${sample}_Split.bam
```

**Input Data:**
- *marked_duplicates.bam (samtools sorted mapping to genome with duplicates marked file, output from step 4a)

For processing human-filtered data
- Homo_sapiens.GRCh38.dna.primary_assembly.fa (genome sequence)

For processing unfiltered data
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa (genome sequence)

**Output Data:**
- *Split.bam (BAM file with reads split at N CIGAR elements and CIGAR strings updated)

---

<br>

## 6. Generate Recalibration Table for Base Quality Score Recalibration

```
gatk --java-options "-Xmx100G" BaseRecalibrator -I $inter_dir/${sample}/${sample}_Split.bam \
  -R $genome_ref/*.fa \
  --known-sites $vcf_ref/Homo_sapiens_assembly38.ens100.known_indels.vcf.gz \
  --known-sites $vcf_ref/dbSNP_v153_ens.vcf.gz \
  -O $inter_dir/${sample}/${sample}_recal_data.table \
  --sequence-dictionary $genome_ref/*.dict
```

**Input Data:**
- *Split.bam (BAM file with reads split at N CIGAR elements and CIGAR strings updated, output from step 5)
- Homo_sapiens_assembly38.ens100.known_indels.vcf.gz (human indels reference with ensembl IDs, output from step 2b)
- dbSNP_v153_ens.vcf.gz (human SNP database with ensembl IDs, output from step 2b)

For processing human-filtered data
- Homo_sapiens.GRCh38.dna.primary_assembly.fa (genome sequence)
- Homo_sapiens.GRCh38.dna.primary_assembly.dict (genome dictionary, output from step 1c)

For processing unfiltered data
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa (genome sequence)
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.dict (genome dictionary, output from step 1c)

**Output Data:**
- *recal_data.table (text file containing the list of arguments used, tables with quantized qualities, tables containing recalibration by read group, quality score, and the optional covariates)

---

<br>

## 7. Evaluate and Compare Base Quality Score Recalibration Tables

```
gatk --java-options "-Xmx100G" AnalyzeCovariates -bqsr $inter_dir/${sample}/${sample}_recal_data.table \
  -plots $inter_dir/${sample}/${sample}_AnalyzeCovariates.pdf
```

**Input Data:**
- *recal_data.table (text file containing the list of arguments used, tables with quantized qualities, tables containing recalibration by read group, quality score, and the optional covariates, output from step 6)

**Output Data:**
- *AnalyzeCovariates.pdf (pdf file containing plots to assess the quality of the recalibration)

---

<br>

## 8. Apply Base Quality Score Recalibration

```
gatk --java-options "-Xmx100G" ApplyBQSR -R $genome_ref/*.fa \
  -I $inter_dir/${sample}/${sample}_Split.bam \
  --bqsr-recal-file $inter_dir/${sample}/${sample}_recal_data.table \
  -O $inter_dir/${sample}/${sample}_BSQR-applied.out.bam
```

**Input Data:**
- *Split.bam (BAM file with reads split at N CIGAR elements and CIGAR strings updated, output from step 5)
- *recal_data.table (text file containing the list of arguments used, tables with quantized qualities, tables containing recalibration by read group, quality score, and the optional covariates, output from step 6)

For processing human-filtered data
- Homo_sapiens.GRCh38.dna.primary_assembly.fa (genome sequence)

For processing unfiltered data
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa (genome sequence)

**Output Data:**
- *BSQR-applied.out.bam (BAM file containing the recalibrated read data)

<br>

---

**In steps 9 - 12 chromosomes are split for parallel processing**

---

<br>

## 9. Call Germline SNPs and Indels via Local Re-assembly of Haplotypes

```
gatk --java-options "-Xmx100G" HaplotypeCaller -R $genome_ref/*.fa \
  -I $inter_dir/${sample}/${sample}_BSQR-applied.out.bam \
  -O $out_dir/${sample}/${chr}/${sample}_${chr}.vcf.gz \
  -ERC GVCF \
  --dont-use-soft-clipped-bases \
  --standard-min-confidence-threshold-for-calling 20 \
  -G StandardAnnotation \
  -G AS_StandardAnnotation \
  --intervals ${chr}
```

**Input Data:**
- *BSQR-applied.out.bam (BAM file containing the recalibrated read data, output from step 8)
- ${chr} (ensembl chromosome identifer)

For processing human-filtered data
- Homo_sapiens.GRCh38.dna.primary_assembly.fa (genome sequence)

For processing unfiltered data
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa (genome sequence)

**Output Data:**
- *.vcf.gz (GVCF file containing raw, unfiltered SNP and indel calls)

---

<br>

## 10. Import VCFs to Genomics Database

```
gatk --java-options "-Xmx40G" GenomicsDBImport -V $out_dir/sample1/${chr}/sample1_${chr}.vcf.gz \
  -V $out_dir/sample2/${chr}/sample2_${chr}.vcf.gz \
  -V $out_dir/sample3/${chr}/sample3_${chr}.vcf.gz \
  <…> \
  -V $out_dir/sample731/${chr}/sample731_${chr}.vcf.gz \
  -V $out_dir/sample732/${chr}/sample732_${chr}.vcf.gz \
  --genomicsdb-workspace-path $out_dir/GVCF_databases/${chr}_database \
  --intervals ${chr}
```

**Input Data:**
- *.vcf.gz (GVCF file containing raw, unfiltered SNP and indel calls, output from step 9)
- ${chr} (ensembl chromosome identifer)

**Output Data:**
- *database (genomics database containing raw, unfiltered SNP and indel calls for all samples that will be joint-genotyped in step 11)

---

<br>

## 11. Perform Joint Genotyping

```
gatk --java-options "-Xmx40G" GenotypeGVCFs -R $genome_ref/*.fa \
  -V gendb://$out_dir/GVCF_databases/${chr}_database \
  -G StandardAnnotation \
  -G AS_StandardAnnotation \
  -O $out_dir/all_samples/${chr}/${chr}_Geno_out.vcf.gz
```

**Input Data:**
- *database (genomics database, output from step 10)
- ${chr} (ensembl chromosome identifer)

For processing human-filtered data
- Homo_sapiens.GRCh38.dna.primary_assembly.fa (genome sequence)

For processing unfiltered data
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa (genome sequence)

**Output Data:**
- *Geno_out.vcf.gz (VCF file containing all jointly genotyped samples)

---

<br>

## 12. Filter Variant Calls Based on Annotations

```
gatk VariantFiltration -R $genome_ref/*.fa \
  -V $out_dir/all_samples/${chr}/${chr}_Geno_out.vcf.gz \
  -O $out_dir/all_samples/${chr}/${chr}_VarFilt_output.vcf.gz \
  --window 35 \
  --cluster 3 \
  --filter-name "FS" \
  --filter "FS > 30.0" \
  --filter-name "QD" \
  --filter "QD < 2.0"
```

**Input Data:**
- *Geno_out.vcf.gz (VCF file containing all jointly genotyped samples, output from step 11)
- ${chr} (ensembl chromosome identifer)

For processing human-filtered data
- Homo_sapiens.GRCh38.dna.primary_assembly.fa (genome sequence)

For processing unfiltered data
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa (genome sequence)

**Output Data:**
- *VarFilt_output.vcf.gz (filtered VCF file containing variants annotated as having passed or failed the indicated `filters`)

---

<br>

## 13. Combine Variant Files From Each Chromosome

```
gatk MergeVcfs --INPUT $out_dir/all_samples/chr1/chr1_VarFilt_output.vcf.gz \
  --INPUT $out_dir/all_samples/chr2/chr2_VarFilt_output.vcf.gz \
  --INPUT $out_dir/all_samples/chr3/chr3_VarFilt_output.vcf.gz \
  <…> \
  --INPUT $out_dir/all_samples/chrX/chrX_VarFilt_output.vcf.gz \
  --INPUT $out_dir/all_samples/chrY/chrY_VarFilt_output.vcf.gz \
  --INPUT $out_dir/all_samples/chrMT/chrMT_VarFilt_output.vcf.gz \
  --OUTPUT $out_dir/all_samples/merged_chr.vcf.gz \
  --SEQUENCE_DICTIONARY $genome_ref/*.dict
```

**Input Data:**
- *VarFilt_output.vcf.gz (filtered VCF file, output from step 12)

For processing human-filtered data
- Homo_sapiens.GRCh38.dna.primary_assembly.dict (genome dictionary, output from step 1c)

For processing unfiltered data
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.dict (genome dictionary, output from step 1c)

**Output Data:**
- *merged_chr.vcf.gz (filtered VCF file containing all chromosomes for all samples sorted according to the `SEQUENCE_DICTIONARY` and by coordinate)

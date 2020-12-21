# GATK Variant Calls Processing Pipeline

> **This page holds an overview and instructions for how COV-IRT generates variant call data from raw Illumina RNA-sequencing data of COVID-19 samples. This pipeline uses the indexed \*Aligned.sortedByCoord_sorted.out.bam files generated from [step 6 of the 'Raw to Aligned Data Pipeline'](Raw_to_Aligned_Data_Pipeline.md#6-sort-and-index-genome-aligned-data). Exact processing commands used for specific cohorts of samples are available in the [Exact_scripts_used](Exact_scripts_used) sub-directory.**  

---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Retrieve and Index Genome Files and Create Reference Dictionary**]()
    - [**1a. Get Genome Files**]()
    - [**1b. Index Genome Files**]()
    - [**1c. Create Reference Dictionaries**]()
  - [**2. Retrieve Variant Call Reference Files and Convert to Match Genome IDs**]()
    - [**2a. Get Variant Call Reference Files**]()
    - [**2b. Convert Variant Call Reference Files to Match Genome IDs then Index**]()
  - [**3. Define Paths to Directories Containing Input, Reference, Intermediate, and Output Files**]()
  - [**4. Mark and Index Duplicate Reads**]()
    - [**4a. Mark Duplicate Reads**]()
    - [**4b. Index Duplicate Reads**]()
  - [**5. Split Reads with N in Cigar**]()
  - [**6. Generate Recalibration Table for Base Quality Score Recalibration**]()
  - [**7. Evaluate and Compare Base Quality Score Recalibration Tables**]()
  - [**8. Apply Base Quality Score Recalibration**]()
  - [**9. Call Germline SNPs and Indels via Local Re-assembly of Haplotypes**]()
  - [**10. Import VCFs to Genomics Database**]()
  - [**11. Perform Joint Genotyping**]()
  - [**12. Filter Variant Calls Based on Annotations**]()
  - [**13. Combine Variant Files**]()

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
wget ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa.gz 
```

Concatenate human and SARS-CoV-2 fasta files - concatenated files used for processing unfiltered data

```
zcat Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa
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

Convert human known indels file

```
java -jar fgbio-1.2.0.jar UpdateVcfContigNames -i Homo_sapiens_assembly38.known_indels.vcf.gz \
  --skip-missing true \
  -d ENS100_MN908947.3.dna.primary_assembly.dict \
  -o Homo_sapiens_assembly38.ens100.known_indels.vcf.gz
```

```
tabix -p vcf Homo_sapiens_assembly38.ens100.known_indels.vcf.gz
```

Convert human SNP database

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
STAR --twopassMode Basic \
  --limitBAMsortRAM 65000000000 \
  --genomeDir /path/to/STAR/genome/directory \
  --outSAMunmapped Within \
  --outFilterType BySJout \
  --outSAMattributes NH HI AS nM NM MD jM jI MC ch \ #same as All
  --outSAMattrRGline ID:flowcell.laneX PL:ILLUMINA PU:flowcell.laneX.${index} LB:${sample} SM:${sample} , ID:flowcell.laneY PL:ILLUMINA PU:flowcell.laneY.${index} LB:${sample} SM:${sample} \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \ # for PE only
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --sjdbScore 1 \
  --readFilesCommand zcat \
  --runThreadN NumberOfThreads \
  --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
  --chimOutJunctionFormat 1 \
  --chimSegmentMin 20 \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM GeneCounts \
  --outSAMheaderHD @HD VN:1.4 SO:coordinate \
  --outFileNamePrefix /path/to/STAR/output/directory/${sample}/${sample}_ \
  --readFilesIn /path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneX_R1.fq.gz,/path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneY_R1.fq.gz /path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneX_R2.fq.gz,/path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneY_R2.fq.gz
```

**Input Data:**
- STAR genome reference (output from step 4b)
- \*flowcell_lane#_R\*.fq.gz (trimmed reads split according to flowcell (i.e. sequencing run) and lane number from step 3)

**Output Data:**
- *Aligned.sortedByCoord.out.bam (sorted mapping to genome)
- *Aligned.toTranscriptome.out.bam (sorted mapping to transcriptome)
- *Chimeric.out.junction (chimerically aligned read data)
- *Chimeric.out.sam (sam file containing chimeric alignments)
- *Log.final.out (log file conting alignment info/stats such as reads mapped, etc)
- *Log.out
- *Log.progress.out
- *ReadsPerGene.out.tab (STAR read counts per gene)
- *SJ.out.tab (high confidence collapsed splice junctions)
- *_STARgenome (directory containing the following:)
  - sjdbInfo.txt
  - sjdbList.out.tab
- *_STARpass1 (directory containing the following:)
  - Log.final.out
  - SJ.out.tab
- *_STARtmp (directory containing the following:)
  - BAMsort (directory containing subdirectories that are empty – this was the location for temp files that were automatically removed after successful completion)

---

<br>

## 9. Call Germline SNPs and Indels via Local Re-assembly of Haplotypes

```
STAR --twopassMode Basic \
  --limitBAMsortRAM 65000000000 \
  --genomeDir /path/to/STAR/genome/directory \
  --outSAMunmapped Within \
  --outFilterType BySJout \
  --outSAMattributes NH HI AS nM NM MD jM jI MC ch \ #same as All
  --outSAMattrRGline ID:flowcell.laneX PL:ILLUMINA PU:flowcell.laneX.${index} LB:${sample} SM:${sample} , ID:flowcell.laneY PL:ILLUMINA PU:flowcell.laneY.${index} LB:${sample} SM:${sample} \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \ # for PE only
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --sjdbScore 1 \
  --readFilesCommand zcat \
  --runThreadN NumberOfThreads \
  --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
  --chimOutJunctionFormat 1 \
  --chimSegmentMin 20 \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM GeneCounts \
  --outSAMheaderHD @HD VN:1.4 SO:coordinate \
  --outFileNamePrefix /path/to/STAR/output/directory/${sample}/${sample}_ \
  --readFilesIn /path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneX_R1.fq.gz,/path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneY_R1.fq.gz /path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneX_R2.fq.gz,/path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneY_R2.fq.gz
```

**Input Data:**
- STAR genome reference (output from step 4b)
- \*flowcell_lane#_R\*.fq.gz (trimmed reads split according to flowcell (i.e. sequencing run) and lane number from step 3)

**Output Data:**
- *Aligned.sortedByCoord.out.bam (sorted mapping to genome)
- *Aligned.toTranscriptome.out.bam (sorted mapping to transcriptome)
- *Chimeric.out.junction (chimerically aligned read data)
- *Chimeric.out.sam (sam file containing chimeric alignments)
- *Log.final.out (log file conting alignment info/stats such as reads mapped, etc)
- *Log.out
- *Log.progress.out
- *ReadsPerGene.out.tab (STAR read counts per gene)
- *SJ.out.tab (high confidence collapsed splice junctions)
- *_STARgenome (directory containing the following:)
  - sjdbInfo.txt
  - sjdbList.out.tab
- *_STARpass1 (directory containing the following:)
  - Log.final.out
  - SJ.out.tab
- *_STARtmp (directory containing the following:)
  - BAMsort (directory containing subdirectories that are empty – this was the location for temp files that were automatically removed after successful completion)

---

<br>

## 10. Import VCFs to Genomics Database

```
STAR --twopassMode Basic \
  --limitBAMsortRAM 65000000000 \
  --genomeDir /path/to/STAR/genome/directory \
  --outSAMunmapped Within \
  --outFilterType BySJout \
  --outSAMattributes NH HI AS nM NM MD jM jI MC ch \ #same as All
  --outSAMattrRGline ID:flowcell.laneX PL:ILLUMINA PU:flowcell.laneX.${index} LB:${sample} SM:${sample} , ID:flowcell.laneY PL:ILLUMINA PU:flowcell.laneY.${index} LB:${sample} SM:${sample} \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \ # for PE only
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --sjdbScore 1 \
  --readFilesCommand zcat \
  --runThreadN NumberOfThreads \
  --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
  --chimOutJunctionFormat 1 \
  --chimSegmentMin 20 \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM GeneCounts \
  --outSAMheaderHD @HD VN:1.4 SO:coordinate \
  --outFileNamePrefix /path/to/STAR/output/directory/${sample}/${sample}_ \
  --readFilesIn /path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneX_R1.fq.gz,/path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneY_R1.fq.gz /path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneX_R2.fq.gz,/path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneY_R2.fq.gz
```

**Input Data:**
- STAR genome reference (output from step 4b)
- \*flowcell_lane#_R\*.fq.gz (trimmed reads split according to flowcell (i.e. sequencing run) and lane number from step 3)

**Output Data:**
- *Aligned.sortedByCoord.out.bam (sorted mapping to genome)
- *Aligned.toTranscriptome.out.bam (sorted mapping to transcriptome)
- *Chimeric.out.junction (chimerically aligned read data)
- *Chimeric.out.sam (sam file containing chimeric alignments)
- *Log.final.out (log file conting alignment info/stats such as reads mapped, etc)
- *Log.out
- *Log.progress.out
- *ReadsPerGene.out.tab (STAR read counts per gene)
- *SJ.out.tab (high confidence collapsed splice junctions)
- *_STARgenome (directory containing the following:)
  - sjdbInfo.txt
  - sjdbList.out.tab
- *_STARpass1 (directory containing the following:)
  - Log.final.out
  - SJ.out.tab
- *_STARtmp (directory containing the following:)
  - BAMsort (directory containing subdirectories that are empty – this was the location for temp files that were automatically removed after successful completion)

---

<br>

## 11. Perform Joint Genotyping

```
STAR --twopassMode Basic \
  --limitBAMsortRAM 65000000000 \
  --genomeDir /path/to/STAR/genome/directory \
  --outSAMunmapped Within \
  --outFilterType BySJout \
  --outSAMattributes NH HI AS nM NM MD jM jI MC ch \ #same as All
  --outSAMattrRGline ID:flowcell.laneX PL:ILLUMINA PU:flowcell.laneX.${index} LB:${sample} SM:${sample} , ID:flowcell.laneY PL:ILLUMINA PU:flowcell.laneY.${index} LB:${sample} SM:${sample} \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \ # for PE only
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --sjdbScore 1 \
  --readFilesCommand zcat \
  --runThreadN NumberOfThreads \
  --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
  --chimOutJunctionFormat 1 \
  --chimSegmentMin 20 \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM GeneCounts \
  --outSAMheaderHD @HD VN:1.4 SO:coordinate \
  --outFileNamePrefix /path/to/STAR/output/directory/${sample}/${sample}_ \
  --readFilesIn /path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneX_R1.fq.gz,/path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneY_R1.fq.gz /path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneX_R2.fq.gz,/path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneY_R2.fq.gz
```

**Input Data:**
- STAR genome reference (output from step 4b)
- \*flowcell_lane#_R\*.fq.gz (trimmed reads split according to flowcell (i.e. sequencing run) and lane number from step 3)

**Output Data:**
- *Aligned.sortedByCoord.out.bam (sorted mapping to genome)
- *Aligned.toTranscriptome.out.bam (sorted mapping to transcriptome)
- *Chimeric.out.junction (chimerically aligned read data)
- *Chimeric.out.sam (sam file containing chimeric alignments)
- *Log.final.out (log file conting alignment info/stats such as reads mapped, etc)
- *Log.out
- *Log.progress.out
- *ReadsPerGene.out.tab (STAR read counts per gene)
- *SJ.out.tab (high confidence collapsed splice junctions)
- *_STARgenome (directory containing the following:)
  - sjdbInfo.txt
  - sjdbList.out.tab
- *_STARpass1 (directory containing the following:)
  - Log.final.out
  - SJ.out.tab
- *_STARtmp (directory containing the following:)
  - BAMsort (directory containing subdirectories that are empty – this was the location for temp files that were automatically removed after successful completion)

---

<br>

## 12. Filter Variant Calls Based on Annotations

```
STAR --twopassMode Basic \
  --limitBAMsortRAM 65000000000 \
  --genomeDir /path/to/STAR/genome/directory \
  --outSAMunmapped Within \
  --outFilterType BySJout \
  --outSAMattributes NH HI AS nM NM MD jM jI MC ch \ #same as All
  --outSAMattrRGline ID:flowcell.laneX PL:ILLUMINA PU:flowcell.laneX.${index} LB:${sample} SM:${sample} , ID:flowcell.laneY PL:ILLUMINA PU:flowcell.laneY.${index} LB:${sample} SM:${sample} \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \ # for PE only
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --sjdbScore 1 \
  --readFilesCommand zcat \
  --runThreadN NumberOfThreads \
  --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
  --chimOutJunctionFormat 1 \
  --chimSegmentMin 20 \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM GeneCounts \
  --outSAMheaderHD @HD VN:1.4 SO:coordinate \
  --outFileNamePrefix /path/to/STAR/output/directory/${sample}/${sample}_ \
  --readFilesIn /path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneX_R1.fq.gz,/path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneY_R1.fq.gz /path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneX_R2.fq.gz,/path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneY_R2.fq.gz
```

**Input Data:**
- STAR genome reference (output from step 4b)
- \*flowcell_lane#_R\*.fq.gz (trimmed reads split according to flowcell (i.e. sequencing run) and lane number from step 3)

**Output Data:**
- *Aligned.sortedByCoord.out.bam (sorted mapping to genome)
- *Aligned.toTranscriptome.out.bam (sorted mapping to transcriptome)
- *Chimeric.out.junction (chimerically aligned read data)
- *Chimeric.out.sam (sam file containing chimeric alignments)
- *Log.final.out (log file conting alignment info/stats such as reads mapped, etc)
- *Log.out
- *Log.progress.out
- *ReadsPerGene.out.tab (STAR read counts per gene)
- *SJ.out.tab (high confidence collapsed splice junctions)
- *_STARgenome (directory containing the following:)
  - sjdbInfo.txt
  - sjdbList.out.tab
- *_STARpass1 (directory containing the following:)
  - Log.final.out
  - SJ.out.tab
- *_STARtmp (directory containing the following:)
  - BAMsort (directory containing subdirectories that are empty – this was the location for temp files that were automatically removed after successful completion)

---

<br>

## 13. Combine Variant Files

```
STAR --twopassMode Basic \
  --limitBAMsortRAM 65000000000 \
  --genomeDir /path/to/STAR/genome/directory \
  --outSAMunmapped Within \
  --outFilterType BySJout \
  --outSAMattributes NH HI AS nM NM MD jM jI MC ch \ #same as All
  --outSAMattrRGline ID:flowcell.laneX PL:ILLUMINA PU:flowcell.laneX.${index} LB:${sample} SM:${sample} , ID:flowcell.laneY PL:ILLUMINA PU:flowcell.laneY.${index} LB:${sample} SM:${sample} \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \ # for PE only
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --sjdbScore 1 \
  --readFilesCommand zcat \
  --runThreadN NumberOfThreads \
  --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
  --chimOutJunctionFormat 1 \
  --chimSegmentMin 20 \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM GeneCounts \
  --outSAMheaderHD @HD VN:1.4 SO:coordinate \
  --outFileNamePrefix /path/to/STAR/output/directory/${sample}/${sample}_ \
  --readFilesIn /path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneX_R1.fq.gz,/path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneY_R1.fq.gz /path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneX_R2.fq.gz,/path/to/split/trimmed/reads/${sample}/${sample}_flowcell_laneY_R2.fq.gz
```

**Input Data:**
- STAR genome reference (output from step 4b)
- \*flowcell_lane#_R\*.fq.gz (trimmed reads split according to flowcell (i.e. sequencing run) and lane number from step 3)

**Output Data:**
- *Aligned.sortedByCoord.out.bam (sorted mapping to genome)
- *Aligned.toTranscriptome.out.bam (sorted mapping to transcriptome)
- *Chimeric.out.junction (chimerically aligned read data)
- *Chimeric.out.sam (sam file containing chimeric alignments)
- *Log.final.out (log file conting alignment info/stats such as reads mapped, etc)
- *Log.out
- *Log.progress.out
- *ReadsPerGene.out.tab (STAR read counts per gene)
- *SJ.out.tab (high confidence collapsed splice junctions)
- *_STARgenome (directory containing the following:)
  - sjdbInfo.txt
  - sjdbList.out.tab
- *_STARpass1 (directory containing the following:)
  - Log.final.out
  - SJ.out.tab
- *_STARtmp (directory containing the following:)
  - BAMsort (directory containing subdirectories that are empty – this was the location for temp files that were automatically removed after successful completion)

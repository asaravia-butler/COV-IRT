# RSEM Counts Processing Pipeline

> **This page holds an overview and instructions for how COV-IRT generates RSEM count data from raw Illumina RNA-sequencing data of COVID-19 samples. This pipeline uses the \*Aligned.toTranscriptome.out.bam files generated from [step 5a of the 'Raw to Aligned Data Pipeline'](Raw_to_Aligned_Data_Pipeline.md#5a-align-reads-to-reference-genome-with-star). Exact processing commands used for specific cohorts of samples are available in the [Exact_scripts_used](Exact_scripts_used) sub-directory.**  

---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Retrieve Genome/Annotation Files and Build RSEM Reference**](#1-retrieve-genomeannotation-files-and-build-rsem-reference)
    - [**1a. Get Genome and Annotation Files**](#1a-get-genome-and-annotation-files)
    - [**1b. Build RSEM Reference**](#1b-build-rsem-reference)
  - [**2. Count Aligned Reads with RSEM and Generate RSEM Counts Tables**](#2-count-aligned-reads-with-rsem-and-generate-rsem-counts-tables)
    - [**2a. Count Aligned Reads with RSEM**](#2a-count-aligned-reads-with-rsem)
    - [**2b. Generate RSEM Gene and Isoform Counts Tables in R**](#2b-generate-rsem-gene-and-isoform-counts-tables-in-r)

---

<br>

# Software used  

|Program|Version*|Relevant Links|
|:------|:------:|:-------------|
|RSEM|`rsem-calculate-expression --version`|[https://github.com/deweylab/RSEM](https://github.com/deweylab/RSEM)|
|Bioconductor|`BiocManager::version()`|[https://bioconductor.org](https://bioconductor.org)|
|tximport|`packageVersion("tximport")`|[https://bioconductor.org/packages/release/bioc/html/tximport.html](https://bioconductor.org/packages/release/bioc/html/tximport.html)|
|tidyverse|`packageVersion("tidyverse")`|[https://www.tidyverse.org](https://www.tidyverse.org)|

>**\*** Exact versions used to process specific cohorts are available in the [Exact_scripts_used](Exact_scripts_used) sub-directory. 

---

<br>

# General processing overview with example commands  

> Exact processing commands used for specific cohorts are provided in the [Exact_scripts_used](Exact_scripts_used) sub-directory.  

---

## 1. Retrieve Genome/Annotation Files and Build RSEM Reference

### 1a. Get Genome and Annotation Files 

> Genome and annotation files used in this pipeline are the same as those used in [step 4 of the 'Raw to Aligned Data Pipeline'](Raw_to_Aligned_Data_Pipeline.md#4-retrieve-genomeannotation-files-and-build-star-reference).

Get human fasta and gtf files from [Ensembl](https://www.ensembl.org/) - used for processing human-filtered data and needed to process unfiltered data

```
wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
```

Get SARS-CoV-2 fasta and gtf files from [Ensembl](https://www.ensembl.org/) - needed to process unfiltered data

```
wget ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa.gz 

wget ftp://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.100.gtf.gz 
```

Concatenate human and SARS-CoV-2 fasta and gtf files - concatenated files used for processing unfiltered data

```
zcat Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa

zcat Homo_sapiens.GRCh38.100.gtf.gz Sars_cov_2.ASM985889v3.100.gtf.gz > Homo_sapiens.GRCh38.100_and_Sars_cov_2.ASM985889v3.100.gtf 
```

<br>

### 1b. Build RSEM Reference  

```
rsem-prepare-reference --gtf /path/to/annotation/gtf/file \
  /path/to/genome/fasta/file \
  /path/to/RSEM/genome/directory/RSEM_ref_prefix

```

**Input Data:** 

For processing human-filtered data
- Homo_sapiens.GRCh38.dna.primary_assembly.fa (genome sequence)
- Homo_sapiens.GRCh38.100.gtf (genome annotation)

For processing unfiltered data
- Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa (genome sequence)
- Homo_sapiens.GRCh38.100_and_Sars_cov_2.ASM985889v3.100.gtf (genome annotation)

**Output Data:**

RSEM genome reference, which consists of the following files:
- RSEM_ref_prefix.chrlist
- RSEM_ref_prefix.grp
- RSEM_ref_prefix.idx.fa
- RSEM_ref_prefix.n2g.idx.fa
- RSEM_ref_prefix.seq
- RSEM_ref_prefix.ti
- RSEM_ref_prefix.transcripts.fa

---

<br>

## 2. Count Aligned Reads with RSEM and Generate RSEM Counts Tables

### 2a. Count Aligned Reads with RSEM

```
rsem-calculate-expression --num-threads NumberOfThreads \
	--alignments \
	--bam \
	--paired-end \
	--seed 12345 \
	--estimate-rspd \
	--no-bam-output \
	--strandedness reverse \
  --append-names \
	/path/to/STAR/output/directory/${sample}/${sample}_Aligned.toTranscriptome.out.bam \
	/path/to/RSEM/genome/directory/RSEM_ref_prefix \
	/path/to/RSEM/counts/output/directory/${sample}

```

**Input Data:**
- RSEM genome reference (output from step 1b)
- *Aligned.toTranscriptome.out.bam (sorted mapping to transcriptome, output from [step 5a of the 'Raw to Aligned Data Pipeline'](Raw_to_Aligned_Data_Pipeline.md#5a-align-reads-to-reference-genome-with-star))

**Output Data:**
- *genes.results (counts per gene)
- *isoforms.results (counts per isoform)
- *stat (directory containing the following stats files)
	- *cnt
	- *model
	- *theta

<br>

### 2b. Generate RSEM Gene and Isoform Counts Tables in R

```R
print("Make RSEM gene and isoform counts tables")
print("")

work_dir="/path/to/directory/containing/samples.txt/file"
counts_dir="/path/to/directory/containing/RSEM/counts/data"

setwd(file.path(work_dir))

### Pull in sample names ###
study <- read.csv(Sys.glob(file.path(work_dir,"samples.txt")), header = FALSE, row.names = 1, stringsAsFactors = TRUE)


##### Import Gene Counts Data
gene_files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)
# reorder the genes.results files to match the ordering of the samples in the samples.txt file
gene_files <- gene_files[sapply(rownames(study), function(x)grep(x, gene_files, value=FALSE, fixed=TRUE))]
names(gene_files) <- rownames(study)
gene.txi.rsem <- tximport(gene_files, type = "rsem", txIn = FALSE, txOut = FALSE)


##### Export unnormalized gene counts table
setwd(file.path(counts_dir))
write.csv(gene.txi.rsem$counts,file='RSEM_Unnormalized_Gene_Counts.csv')


##### Import Isoform Counts Data
isoform_files <- list.files(file.path(counts_dir),pattern = ".isoforms.results", full.names = TRUE)
# reorder the isoforms.results files to match the ordering of the samples in the samples.txt file
isoform_files <- isoform_files[sapply(rownames(study), function(x)grep(x, isoform_files, value=FALSE, fixed=TRUE))]
names(isoform_files) <- rownames(study)
isoform.txi.rsem <- tximport(isoform_files, type = "rsem", txIn = FALSE, txOut = FALSE)


##### Export unnormalized isoform counts table
setwd(file.path(counts_dir))
write.csv(isoform.txi.rsem$counts,file='RSEM_Unnormalized_Isoform_Counts.csv')


## print session info ##
print("Session Info below: ")
print("")
sessionInfo()

```

**Input Data:**
- *genes.results (RSEM gene counts from step 2a)
- *isoforms.results (RSEM isoform counts from step 2a)
- samples.txt (text file containing a single column list of all samples)

**Output Data:**
- RSEM_Unnormalized_Gene_Counts.csv (table containing RSEM gene counts for all samples)
- RSEM_Unnormalized_Isoform_Counts.csv (table containing RSEM isoform counts for all samples)

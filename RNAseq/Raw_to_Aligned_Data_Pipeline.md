# Raw to Aligned Data Processing Pipeline

> **This page holds an overview and instructions for how COV-IRT generates alignment and STAR count data from raw Illumina RNA-sequencing data of COVID-19 samples. Exact processing commands used for specific cohorts of samples are available in the [Exact_scripts_used](Exact_scripts_used) sub-directory.**  

---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Raw Data QC**]()
    - [**1a. Raw Data QC**](#1a-raw-data-qc)
    - [**1b. Compile Raw Data QC**](#1b-compile-raw-data-qc)
  - [**2. Trim/Filter Raw Data and Trimmed Data QC**]()
    - [**2a. Trim/Filter Raw Data**](#2a-trimfilter-raw-data)
    - [**2b. Trimmed Data QC**](#2b-trimmed-data-qc)
    - [**2c. Compile Trimmed Data QC**](#2c-compile-trimmed-data-qc)
  - [**3. Split Fastq Files Based on Sequencing Run/Lane**]()
  - [**4. Retrieve Genome/Annotation Files and Build STAR Reference**]()
    - [**4a. Get Genome and Annotation Files**]()
    - [**4b. Build STAR Reference**]()
  - [**5. Align Reads to Reference Genome with STAR**](#5-align-reads-to-reference-genome-with-star)
  - [**6. Sort and Index Genome-Aligned Data**]()
    - [**6a. Sort Genome-Aligned Data**]()
    - [**6b. Index Sorted Genome-Aligned Data**]()

---

<br>

# Software used  

|Program|Version*|Relevant Links|
|:------|:------:|:-------------|
|FastQC|`fastqc -v`|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|`multiqc -v`|[https://multiqc.info/](https://multiqc.info/)|
|Trimmomatic|`trimmomatic -version`|[http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)|
|gdc-fastq-splitter|`gdc-fastq-splitter --version`|[https://github.com/kmhernan/gdc-fastq-splitter](https://github.com/kmhernan/gdc-fastq-splitter)|
|STAR|`STAR --version`|[https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)|
|Samtools|`samtools --version`|[http://www.htslib.org/](http://www.htslib.org/)|

>**\*** Exact versions used to process specific cohorts are available in the [Exact_scripts_used](Exact_scripts_used) sub-directory. 

---

<br>

# General processing overview with example commands  

> Exact processing commands used for specific cohorts are provided in the [Exact_scripts_used](Exact_scripts_used) sub-directory.  

---

## 1. Raw Data QC

### 1a. Raw Data QC  

```
fastqc -o /path/to/raw_fastqc/output/directory *.fastq.gz
```

**Input Data:**
- *fastq.gz (raw reads)

**Output Data:**
- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

<br>

### 1b. Compile Raw Data QC  

```
multiqc -n raw_multiqc -o /path/to/raw_multiqc/output/directory /path/to/directory/containing/raw_fastqc/files
```

**Input Data:**
- *fastqc.zip (FastQC data from step 1a)

**Output Data:**
- raw_multiqc.html (multiqc report)
- raw_multiqc_data (directory containing multiqc data)

<br>

---

## 2. Trim/Filter Raw Data and Trimmed Data QC

### 2a. Trim/Filter Raw Data  

```
trimmomatic PE \
  -threads NumberOfThreads \
  -phred33 \
  -trimlog /path/to/trimming/log/outputs/${sample}_trimming.log \
  -summary /path/to/trimming/log/outputs/${sample}_trimming_summary.txt \
  -validatePairs \
  /path/to/raw/reads/${sample}.R1.fastq.gz \ 
  /path/to/raw/reads/${sample}.R2.fastq.gz \
  /path/to/ouput/trimmed/reads/${sample}_R1_P_trimmed.fq.gz \
  /path/to/ouput/trimmed/reads/${sample}_R1_U_trimmed.fq.gz \
  /path/to/ouput/trimmed/reads/${sample}_R2_P_trimmed.fq.gz \
  /path/to/ouput/trimmed/reads/${sample}_R2_U_trimmed.fq.gz \
  ILLUMINACLIP:/path/to/trimmomatic/adapter/fasta/files/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
  LEADING:20 \
  TRAILING:20 \
  MINLEN:15
```

**Input Data:**
- *fastq.gz (raw reads)

**Output Data:**
- *_P_trimmed.fq.gz (paired trimmed reads)
- *_U_trimmed.fq.gz (unpaired trimmed reads)
- *trimming.log (trimming log file)
- *trimming_summary.txt (trimming summary file)

<br>

### 2b. Trimmed Data QC  

```
fastqc -o /path/to/trimmed_fastqc/output/directory *trimmed.fq.gz
```

**Input Data:**
- *_P_trimmed.fq.gz (paired trimmed reads from step 2a)
- *_U_trimmed.fq.gz (unpaired trimmed reads from step 2a)

**Output Data:**
- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

<br>

### 2c. Compile Trimmed Data QC  

```
multiqc -n trimmed_multiqc -o /path/to/trimmed_multiqc/output/directory /path/to/directory/containing/trimmed_fastqc/files
```

**Input Data:**
- *fastqc.zip (FastQC data from step 2b)

**Output Data:**
- trimmed_multiqc.html (multiqc report)
- trimmed_multiqc_data (directory containing multiqc data)

<br>

---

## 3. Split Fastq Files Based on Sequencing Run/Lane


<br>

---

## 4. Build STAR Reference  

```
STAR --runThreadN NumberOfThreads \
  --runMode genomeGenerate \
  --limitGenomeGenerateRAM 55000000000 \
  --genomeSAindexNbases 14 \
  --genomeDir /path/to/STAR/genome/directory \
  --genomeFastaFiles /path/to/genome/fasta/file \
  --sjdbGTFfile /path/to/annotation/gtf/file \
  --sjdbOverhang ReadLength-1

```

**Parameter Definitions:**

* `--runThreadN` – number of threads available on server node to create STAR reference
* `--runMode` - instructs STAR to run genome indices generation job
* `--limitGenomeGenerateRAM` - maximum RAM available (in bytes) to generate STAR reference, at least 35GB are needed for mouse and the example above shows 55GB
* `--genomeSAindexNbases` - length (in bases) of the SA pre-indexing string, usually between 10 and 15. Longer strings require more memory but allow for faster searches. This value should be scaled down for smaller genomes (like bacteria) to min(14, log2(GenomeLength)/2 - 1). For example, for a 1 megaBase genome this value would be 9.
* `--genomeDir` - specifies the path to the directory where the STAR reference will be stored. At least 100GB of available disk space is required for mammalian genomes.
* `--genomeFastaFiles` - specifies one or more fasta file(s) containg the genome reference sequences
* `--sjdbGTFfile` – specifies the file(s) containg annotated transcripts in the standard gtf format
* `--sjdbOverhang` - indicates the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. The length should be one less than the length of the reads.

**Input Data:**
- *.fasta (genome sequence\#)
- *.gtf (genome annotation\#)

\#See document(s) in the [GeneLab_Reference_and_Annotation_Files](GeneLab_Reference_and_Annotation_Files) sub-directory for a list of the ensembl fasta genome sequences and associated gtf annotation files used to generate the RNAseq processed data available in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects).

**Output Data:**

STAR genome reference, which consists of the following files:
- chrLength.txt
- chrNameLength.txt
- chrName.txt
- chrStart.txt
- exonGeTrInfo.tab
- exonInfo.tab
- geneInfo.tab
- Genome
- genomeParameters.txt
- SA
- SAindex
- sjdbInfo.txt
- sjdbList.fromGTF.out.tab
- sjdbList.out.tab
- transcriptInfo.tab

<br>

---

## 4. Align reads to reference genome with STAR

```
STAR --twopassMode Basic \
	--limitBAMsortRAM 65000000000 \
	--genomeDir /path/to/STAR/genome/directory \
	--outSAMunmapped Within \
	--outFilterType BySJout \
	--outSAMattributes NH HI AS NM MD MC \
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
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM \
	--outSAMheaderHD @HD VN:1.4 SO:coordinate \
	--outFileNamePrefix /path/to/STAR/output/directory/<sample_id> \
	--readFilesIn /path/to/trimmed_forward_reads \
	/path/to/trimmed_reverse_reads # only needed for PE studies

```

**Parameter Definitions:**

* `--twopassMode` – specifies 2-pass mapping mode; the `Basic` option instructs STAR to perform the 1st pass mapping, then automatically extract junctions, insert them into the genome index, and re-map all reads in the 2nd mapping pass
* `--limitBAMsortRAM` - maximum RAM available (in bytes) to sort the bam files, the example above indicates 65GB
* `--genomeDir` - specifies the path to the directory where the STAR reference is stored
* `--outSAMunmapped` - specifies ouput of unmapped reads in the sam format; the `Within` option instructs STAR to output the unmapped reads within the main sam file
* `--outFilterType` - specifies the type of filtering; the `BySJout` option instructs STAR to keep only those reads that contain junctions that passed filtering in the SJ.out.tab output file
* `--outSAMattributes` - list of desired sam attributes in the order desired for the output sam file; sam attribute descriptions can be found [here](https://samtools.github.io/hts-specs/SAMtags.pdf)
* `--outFilterMultimapNmax` – specifies the maximum number of loci the read is allowed to map to; all alignments will be output only if the read maps to no more loci than this value
* `--outFilterMismatchNmax` - maximum number of mismatches allowed to be included in the alignment output
* `--outFilterMismatchNoverReadLmax` - ratio of mismatches to read length allowed to be included in the alignment output; the `0.04` value indicates that up to 4 mismatches are allowed per 100 bases
* `--alignIntronMin` - minimum intron size; a genomic gap is considered an intron if its length is equal to or greater than this value, otherwise it is considered a deletion
* `--alignIntronMax` - maximum intron size
* `--alignMatesGapMax` - maximum genomic distance (in bases) between two mates of paired-end reads; this option should be removed for single-end reads
* `--alignSJoverhangMin` - minimum overhang (i.e. block size) for unannotated spliced alignments
* `--alignSJDBoverhangMin` - minimum overhang (i.e. block size) for annotated spliced alignments
* `--sjdbScore` - additional alignment score for alignments that cross database junctions
* `--readFilesCommand` - specifies command needed to interpret input files; the `zcat` option indicates input files are compressed with gzip and zcat will be used to uncompress the gzipped input files
* `--runThreadN` - indicates the number of threads to be used for STAR alignment and should be set to the number of available cores on the server node
* `--outSAMtype` - specifies desired output format; the `BAM SortedByCoordinate` options specify that the output file will be sorted by coordinate and be in the bam format
* `--quantMode` - specifies the type(s) of quantification desired; the `TranscriptomeSAM` option instructs STAR to output a separate sam/bam file containing alignments to the transcriptome
* `--outSAMheaderHD` - indicates a header line for the sam/bam file
* `--outFileNamePrefix` - specifies the path to and prefix for the output file names; for GeneLab the prefix is the sample id
* `--readFilesIn` - path to input read 1 (forward read) and read 2 (reverse read); for paired-end reads, read 1 and read 2 should be separated by a space; for single-end reads only read 1 should be indicated 


**Input Data:**
- STAR genome reference (output from Step 3)
- *fastq.gz (trimmed reads)

**Output Data:**
- *Aligned.sortedByCoord.out.bam# (sorted mapping to genome)
- *Aligned.toTranscriptome.out.bam# (sorted mapping to transcriptome)
- *Log.final.out# (log file conting alignment info/stats such as reads mapped, etc)
- *Log.out
- *Log.progress.out
- *SJ.out.tab\#
- *_STARgenome (directory containing the following:)
  - sjdbInfo.txt
  - sjdbList.out.tab
- *_STARpass1 (directory containing the following:)
  - Log.final.out
  - SJ.out.tab
- *_STARtmp (directory containing the following:)
  - BAMsort (directory containing subdirectories that are empty – this was the location for temp files that were automatically removed after successful completion)

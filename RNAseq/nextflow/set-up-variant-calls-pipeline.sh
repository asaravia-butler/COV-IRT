#!/usr/bin/env bash

# TODO: Clean-up and document

set -ex

REPO_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd | xargs dirname | xargs dirname )"
DATA_DIR="$( dirname $REPO_DIR )/COV-IRT-Data"

mkdir -p $DATA_DIR
pushd $DATA_DIR

if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" ]; then
    wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
fi
if [ ! -f "Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz" ]; then
    wget ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz
fi
if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa" ]; then
    zcat Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz \
	 > Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa
fi
if [ ! -f "Homo_sapiens_assembly38.known_indels.vcf.gz" ]; then
    wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
fi
if [ ! -f "Homo_sapiens_assembly38.known_indels.vcf.gz.tbi" ]; then
    wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
fi
if [ ! -f "00-All.vcf.gz" ]; then
    wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
fi
if [ ! -f "00-All.vcf.gz.tbi" ]; then
    wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi
fi

if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
    gunzip -k Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
fi

conda activate COVIRT_GATK

if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai" ]; then
    samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
fi

if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa.fai" ]; then
    samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa
fi

if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.dict" ]; then
    gatk --java-options "-Xmx100G" CreateSequenceDictionary \
	 -R Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	 -O Homo_sapiens.GRCh38.dna.primary_assembly.dict
fi

if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.dict" ]; then
    gatk --java-options "-Xmx100G" CreateSequenceDictionary \
	 -R Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa \
	 -O Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.dict
fi

conda deactivate

conda activate COVIRT_fgbio

if [ ! -f "Homo_sapiens_assembly38.ens100.known_indels.vcf.gz" ]; then
    fgbio UpdateVcfContigNames -i Homo_sapiens_assembly38.known_indels.vcf.gz \
	  --skip-missing=true \
	  -d $REPO_DIR/RNAseq/Reference_Files/ENS100_MN908947.3.dna.primary_assembly.dict \
	  -o Homo_sapiens_assembly38.ens100.known_indels.vcf.gz
fi

if [ ! -f "dbSNP_v153_ens.vcf.gz" ]; then
    fgbio UpdateVcfContigNames -i 00-All.vcf.gz \
	  --skip-missing=true \
	  -d $REPO_DIR/RNAseq/Reference_Files/ENS100_MN908947.3.dna.primary_assembly.dict \
	  -o dbSNP_v153_ens.vcf.gz
fi

conda deactivate

popd

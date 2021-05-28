#!/usr/bin/env bash

# TODO: Clean-up and document

set -ex

REPO_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd | xargs dirname | xargs dirname )"
# DATA_DIR="$( dirname $REPO_DIR )/COV-IRT-Data"
DATA_DIR="/nobackupp13/rleclai2/COV-IRT-Data"

mkdir -p $DATA_DIR
pushd $DATA_DIR

if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" ]; then
    wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
fi
if [ ! -f "Homo_sapiens.GRCh38.100.gtf.gz" ]; then
    wget ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
fi
if [ ! -f "Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz" ]; then
    wget ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz
fi
if [ ! -f "Sars_cov_2.ASM985889v3.101.gtf.gz" ]; then
    wget ftp://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz
fi

if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa" ]; then
    zcat Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz \
	 > Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa
fi
if [ ! -f "Homo_sapiens.GRCh38.100_and_Sars_cov_2.ASM985889v3.100.gtf" ]; then
    zcat Homo_sapiens.GRCh38.100.gtf.gz Sars_cov_2.ASM985889v3.101.gtf.gz \
	 > Homo_sapiens.GRCh38.100_and_Sars_cov_2.ASM985889v3.100.gtf 
fi

if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
    gunzip -k Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
fi
if [ ! -f "Homo_sapiens.GRCh38.100.gtf" ]; then
    gunzip -k Homo_sapiens.GRCh38.100.gtf.gz
fi

conda activate COVIRT_fastq_to_alignment

NumberOfThreads=16
# TODO: Check read length
ReadLengthM1=150

if [ ! -d "filtered" ]; then
    STAR --runThreadN ${NumberOfThreads} \
	 --runMode genomeGenerate \
	 --limitGenomeGenerateRAM 55000000000 \
	 --genomeSAindexNbases 14 \
	 --genomeDir filtered \
	 --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	 --sjdbGTFfile Homo_sapiens.GRCh38.100.gtf \
	 --sjdbOverhang ${ReadLengthM1}
fi

if [ ! -d "unfiltered" ]; then
    STAR --runThreadN ${NumberOfThreads} \
	 --runMode genomeGenerate \
	 --limitGenomeGenerateRAM 55000000000 \
	 --genomeSAindexNbases 14 \
	 --genomeDir unfiltered \
	 --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa \
	 --sjdbGTFfile Homo_sapiens.GRCh38.100_and_Sars_cov_2.ASM985889v3.100.gtf \
	 --sjdbOverhang ${ReadLengthM1}
fi

conda deactivate

popd

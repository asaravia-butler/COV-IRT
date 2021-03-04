#!/usr/bin/env bash

# TODO: Clean-up and document

NumberOfThreads=64
# TODO: Check read length
ReadLengthM1=150

STAR --runThreadN ${NumberOfThreads} \
     --runMode genomeGenerate \
     --limitGenomeGenerateRAM 55000000000 \
     --genomeSAindexNbases 14 \
     --genomeDir filtered \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile Homo_sapiens.GRCh38.100.gtf \
     --sjdbOverhang ${ReadLengthM1}

STAR --runThreadN ${NumberOfThreads} \
     --runMode genomeGenerate \
     --limitGenomeGenerateRAM 55000000000 \
     --genomeSAindexNbases 14 \
     --genomeDir unfiltered \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa \
     --sjdbGTFfile Homo_sapiens.GRCh38.100_and_Sars_cov_2.ASM985889v3.100.gtf \
     --sjdbOverhang ${ReadLengthM1}

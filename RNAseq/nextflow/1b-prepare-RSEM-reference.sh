#!/usr/bin/env sh

# rsem-prepare-reference doesn't seem to be okay with gzipped file

rsem-prepare-reference --gtf Homo_sapiens.GRCh38.100_and_Sars_cov_2.ASM985889v3.100.gtf Homo_sapiens.GRCh38.dna.primary_assembly_and_Sars_cov_2.ASM985889v3.dna.primary_assembly.MN908947.3.fa Homo_sapiens_and_Sars_cov_2_RSEMref

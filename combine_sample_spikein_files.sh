#!/bin/bash

# Input:
sample_vcf=$1
spikein_vcf=$2
combined_vcf_gz=$3
sample_gtf=$4
spikein_gtf=$5
combined_gtf=$6
sample_fasta_a1=$7
spikein_fasta_a1=$8
combined_fasta_a1=$9
sample_fasta_a2=${10}
spikein_fasta_a2=${11}
combined_fasta_a2=${12}

# Combine VCF files:
bcftools concat ${sample_vcf} ${spikein_vcf} -O v | bcftools sort -O v | bgzip > ${combined_vcf_gz}
tabix -p vcf ${combined_vcf_gz}

# Combine GTF files:
cat ${sample_gtf} ${spikein_gtf} | sort -k1,1 -k4,4n > ${combined_gtf}

# Combine FASTA files:
cat ${sample_fasta_a1} ${spikein_fasta_a1} > ${combined_fasta_a1}
samtools faidx ${combined_fasta_a1}
cat ${sample_fasta_a2} ${spikein_fasta_a2} > ${combined_fasta_a2}
samtools faidx ${combined_fasta_a2}

# Message:
echo "Combination completed:"
echo "  Combined VCF: ${combined_vcf_gz}"
echo "  Combined GTF: ${combined_gtf}"
echo "  Combined FASTA A1: ${combined_fasta_a1}"
echo "  Combined FASTA A2: ${combined_fasta_a2}"

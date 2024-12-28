#!/bin/bash

# Input:
suffix=$1
fastaA1=$2
fastaA2=$3
vcf=$4
gtf=$5
output_dir=$6

# New filenames and directory:
names_txt="${output_dir}/old_and_new_names.txt"
renamed_fastaA1="${output_dir}/$(basename ${fastaA1%.fa}).chr_renamed.fa"
renamed_fastaA2="${output_dir}/$(basename ${fastaA2%.fa}).chr_renamed.fa"
renamed_vcf="${output_dir}/$(basename ${vcf%.vcf}).chr_renamed.vcf.gz"
renamed_gtf="${output_dir}/$(basename ${gtf%.gtf}).chr_renamed.gtf"
mkdir -p "$output_dir"

# Renaming VCF: 
awk -v suffix="$suffix" '{print $1 "\t" $1 suffix}' ${fastaA1}.fai > ${names_txt}
bcftools annotate --rename-chrs ${names_txt} ${vcf} -Oz | \
bcftools annotate --rename-chrs ${names_txt} -Oz | \
bcftools sort -Oz -o ${renamed_vcf} 
tabix ${renamed_vcf}

# Renaming FASTA:
sed "/^>/ s/^\(>\S*\)/\1${suffix}/" ${fastaA1} > ${renamed_fastaA1}
sed "/^>/ s/^\(>\S*\)/\1${suffix}/" ${fastaA2} > ${renamed_fastaA2}

# Renaming GTF:
sed 's/ /xlifehack/g' ${gtf} | \
  awk -v suffix="$suffix" 'BEGIN {OFS="\t"}; !/^#/ {$1 = $1 suffix; print $0;}' | \
  sed 's/xlifehack/ /g' > ${renamed_gtf}

# Message:
echo "Chromosomes renaming completed. Results:"
echo "  (suffix = ${suffix})"
echo "  FASTA A1: ${renamed_fastaA1}"
echo "  FASTA A2: ${renamed_fastaA2}"
echo "  VCF: ${renamed_vcf}"
echo "  GTF: ${renamed_gtf}"


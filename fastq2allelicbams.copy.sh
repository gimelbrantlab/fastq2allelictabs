#!/bin/bash
#SBATCH -p pool
#SBATCH -c 8
#SBATCH --mem-per-cpu=20G
#SBATCH -o /home/amendelevich/logs/palladium/13apr23/fastq2albam_%A_%a.out
#SBATCH -e /home/amendelevich/logs/palladium/13apr23/fastq2albam_%A_%a.err

module load gcc/5.3.0 
module load samtools/1.12 bcftools/1.7 htslib/1.7 python/3.6.10

nthreds=8
d_f2t=/home/amendelevich/tools/fastq2allelictabs/

# Table to fill:
infoTab=$1

# ----------
# Take info:
i=${SLURM_ARRAY_TASK_ID}
#-# --array=n-m 

infoVect=(`grep -w "^"$i $infoTab`)
# 13 mandatory fields (0-12), so:
sampleID=${infoVect[1]}
refSample=${infoVect[4]}
refSpikeIn=${infoVect[5]}
chimericName=${infoVect[7]}
refDir1=${infoVect[8]}
refDir2=${infoVect[9]}
fastq1Path=${infoVect[10]}
fastq2Path=${infoVect[11]}
d=${infoVect[12]}

refDirs=( 0 $refDir1 $refDir2 ) 

# -----------
# Processing:
alD=$d/alignments/$sampleID/
mkdir -p $alD

for A in 1 2; do
  vcf=$(ls ${refDirs[$A]}/*vcf)

  alprefix=$alD/$sampleID".on_"$chimericName".on_Allele"$A"."
  albam=$alprefix"Aligned.sortedByCoord.out.bam"
  fsam=${albam::-3}"filtered_vW1.sam"
  fsamNSort=${fsam::-3}"Nsort.sam"

  # STAR alignment:
  echo "Aligning $A ..."

  /home/amendelevich/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR \
      --readFilesIn $fastq1Path $fastq2Path \
      --readFilesCommand gunzip -c \
      --outFileNamePrefix $alprefix \
      --runThreadN $nthreds \
      --outSAMtype BAM SortedByCoordinate \
      --varVCFfile $vcf \
      --outSAMattributes NH HI AS nM vA vG \
      --genomeDir ${refDirs[$A]} \
      --outFilterMultimapNmax 1 \
      --waspOutputMode SAMtag \
      --outSAMattrRGline ID:Allele$A
      #--quantMode GeneCounts

  # Allelic reads filtering: 
  samtools view --no-PG -H $albam  > $fsam
  samtools view --no-PG $albam | grep "vW:i:1" >> $fsam
  samtools sort --no-PG -n -o $fsamNSort $fsam
  rm $fsam
done

# Allele separation:
echo "Separating alleles..."

fsamNSort1=$alD/$sampleID".on_"$chimericName".on_Allele1.Aligned.sortedByCoord.out.filtered_vW1.Nsort.sam"
fsamNSort2=$alD/$sampleID".on_"$chimericName".on_Allele2.Aligned.sortedByCoord.out.filtered_vW1.Nsort.sam"

python3 $d_f2t/fastq_to_allelic_counts_tabs/alleleseparation.py \
        --samA1 $fsamNSort1 --samA2 $fsamNSort2 \
        --obase $sampleID"_"$chimericName --odir $alD \
        --allele1 "Allele1" --allele2 "Allele2" \
        --paired 1 

for sam in $(ls $alD/*.sam); do
  bam=${sam::-3}bam
  samtools view --no-PG -o $bam $sam
done

rm $alD/*.sam

# -------------
# SNP coverage:

# #############
# #############

# ----------------
# chrom cov stats:

fbamRef=$alD/$sampleID"_"$chimericName".Allele1.bam"
fbamAlt=$alD/$sampleID"_"$chimericName".Allele2.bam"
albam1=$alD/$sampleID".on_"$chimericName".on_Allele1.Aligned.sortedByCoord.out.bam"
albam2=$alD/$sampleID".on_"$chimericName".on_Allele2.Aligned.sortedByCoord.out.bam"

chr_counts=$sampleID"_"$chimericName".chr_counts"
bams_to_chr_count=( $albam1 $albam2 $fbamRef $fbamAlt )
suffs=(".aligned_on_Allele1.tsv" ".aligned_on_Allele2.tsv" ".Allele1.tsv" ".Allele2.tsv")
for ib in {0..3}; do
    samtools view --no-PG ${bams_to_chr_count[$ib]} | cut -f3 | sort -V | uniq -c | sed -e 's/^ *//' -e 's/ /\t/' > $alD/$chr_counts${suffs[$ib]}
done


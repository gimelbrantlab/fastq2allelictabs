#!/bin/bash
#SBATCH -p pool
#SBATCH -c 8
#SBATCH --mem-per-cpu=20G
#SBATCH -o /home/amendelevich/logs/palladium/13apr23/fastq2albam_%A_%a.out
#SBATCH -e /home/amendelevich/logs/palladium/13apr23/fastq2albam_%A_%a.err

module load samtools/1.12 bcftools/1.7 htslib/1.7 python/3.6.10 GATK/4.0.1.0

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
  rsam=${albam::-3}"nonallelic_wasp.sam"
  rsamNSort=${rsam::-3}"Nsort.sam"

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
  samtools view --no-PG -H $albam  > $rsam
  samtools view --no-PG $albam | grep -v "vW:i:1" >> $rsam
  samtools sort --no-PG -n -o $rsamNSort $rsam
  rm $rsam
done

  # who's done here?
  fsamNSort1=$alD/$sampleID".on_"$chimericName".on_Allele1.Aligned.sortedByCoord.out.filtered_vW1.Nsort.sam"
  fsamNSort2=$alD/$sampleID".on_"$chimericName".on_Allele2.Aligned.sortedByCoord.out.filtered_vW1.Nsort.sam"
  rsamNSort1=$alD/$sampleID".on_"$chimericName".on_Allele1.Aligned.sortedByCoord.out.nonallelic_wasp.Nsort.sam"
  rsamNSort2=$alD/$sampleID".on_"$chimericName".on_Allele2.Aligned.sortedByCoord.out.nonallelic_wasp.Nsort.sam"

# Allele separation for BAMs with allelic reads:
echo "Separating alleles..."

python3 $d_f2t/fastq_to_allelic_counts_tabs/alleleseparation.py \
        --samA1 $fsamNSort1 --samA2 $fsamNSort2 \
        --obase $sampleID"_"$chimericName --odir $alD \
        --allele1 "Allele1" --allele2 "Allele2" \
        --paired 1 

  # who's done here?
  asHapNAbam=$alD/$sampleID"_"$chimericName".ambiguous.bam"
  asHap1bam=$alD/$sampleID"_"$chimericName".Allele1.bam"
  asHap2bam=$alD/$sampleID"_"$chimericName".Allele2.bam"

#  ... cleaning and compressing:
for sam in $(ls $alD/*.sam); do
  bam=${sam::-3}bam
  samtools view --no-PG -o $bam $sam
done

rm $alD/*.sam


# Merging BAMs with non-allelic reads:
echo "Merging non-allelic bams..."

rsamNSort=$sampleID"_"$chimericName".nonallelic_wasp.sam"
rbamNSort=$sampleID"_"$chimericName".nonallelic_wasp.bam"
bamHapNA=$sampleID"_"$chimericName".nonallelic_all.bam"

#   -- merge non-allelic bams
### python3 $d_f2t/fastq_to_allelic_counts_tabs/merge_intersecting_sams.py 
### ${rsamNSort1} ${rsamNSort2}
############################################################
samtools view --no-PG -o $rbamNSort $rsamNSort
rm $rsamNSort
#   -- merge non-alelic and ambigous
samtools merge -n -o ${bamHapNA}  ${rbamNSort} ${asHapNAbam} 
#   -- change RG to HapNA in any read 
################ RG CHANGING X1 ############################

# Join assigned Allele1, Allele2 and Unassigned, and index them
#   -- change RG to Hap1, Hap2 in assigned to Allele1 and Allele2 bams
tmpH1bam=${asHap1bam}".tmp.bam"
tmpH2bam=${asHap2bam}".tmp.bam"
gatk AddOrReplaceReadGroups \
       -I $asHap1bam -O $tmpH1bam --RGID Hap1 \
       --RGLB xx --RGPL ILLUMINA --RGPU xx  --RGSM xx
# WRONG: https://www.biostars.org/p/249376/
################ RG CHANGING X2 ############################
#   -- merge with unassigned
IGVbam=$alD/$sampleID"_"$chimericName".IGV.bam"
samtools merge -n -o ${IGVbam}  ${bamHapNA} ${tmpH1bam} ${tmpH2bam}
rm ${tmpH1bam}; rm ${tmpH2bam}
#   -- index resulting bam
samtools index ${IGVbam}


# ----------------
# SNP coverage cmp 
# (to be added?)

# #############
# #############

# ----------------
# --- preSTATS ---
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


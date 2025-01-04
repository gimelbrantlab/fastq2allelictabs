#!/bin/bash

# EXAMPLE:  
#   sbatch $scripts_d/fastq2allelicbams.sh $TAB $i_col1
#   sbatch ~/my_dir/fastq2allelicbams.sh /path/to/table_info.tsv 1

# Table to fill:
infoTab=$1
# Value from the first column:
i=$2

nthreds=4 # subject to change if needed.

# ----------
# Take info:

infoVect=($(awk -v i=$i '{if ($1==i) {print}}' $infoTab))
# 15 mandatory fields (0-14), so:
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

  STAR \
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

# ^^ what's created:
fsamNSort1=$alD/$sampleID".on_"$chimericName".on_Allele1.Aligned.sortedByCoord.out.filtered_vW1.Nsort.sam"
fsamNSort2=$alD/$sampleID".on_"$chimericName".on_Allele2.Aligned.sortedByCoord.out.filtered_vW1.Nsort.sam"
rsamNSort1=$alD/$sampleID".on_"$chimericName".on_Allele1.Aligned.sortedByCoord.out.nonallelic_wasp.Nsort.sam"
rsamNSort2=$alD/$sampleID".on_"$chimericName".on_Allele2.Aligned.sortedByCoord.out.nonallelic_wasp.Nsort.sam"

# Allele separation for BAMs with allelic reads:
echo "Separating alleles..."

python3 /opt/fastq2allelictabs/fastq_to_allelic_counts_tabs/alleleseparation.py \
        --samA1 $fsamNSort1 --samA2 $fsamNSort2 \
        --obase $sampleID"_"$chimericName --odir $alD \
        --allele1 "Allele1" --allele2 "Allele2" \
        --paired 1 

# ^^ what's created:
asHapNAbam=$alD/$sampleID"_"$chimericName".ambiguous.bam"
asHap1bam=$alD/$sampleID"_"$chimericName".Allele1.bam"
asHap2bam=$alD/$sampleID"_"$chimericName".Allele2.bam"

#  ... cleaning and compressing:
for sam in $(ls $alD/*.sam); do
  bam=${sam::-3}bam
  samtools view --no-PG -o $bam $sam
done

rm $alD/*.sam


# ---------------------
# chrom coverage stats:

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





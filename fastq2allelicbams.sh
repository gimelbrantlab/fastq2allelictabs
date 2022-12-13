#!/bin/bash
#SBATCH -p pool
#SBATCH -c 4
#SBATCH --mem-per-cpu=16G
#SBATCH -o /home/amendelevich/logs/bulkSpike/Sep22/fastq2albam_%A_%a.out
#SBATCH -e /home/amendelevich/logs/bulkSpike/Sep22/fastq2albam_%A_%a.err

module load gcc/5.3.0 
module load samtools/1.12 bcftools/1.7 htslib/1.7 python/3.6.10

nthreds=4

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
refDir=${infoVect[9]}
fastq1Path=${infoVect[10]}
fastq2Path=${infoVect[11]}
d=${infoVect[12]}
 
vcf=$(ls $refDir/*vcf)

# -----------
# Processing:
alD=$d/alignments/$sampleID/
mkdir -p $alD
alprefix=$alD/$sampleID".on_"$refSample"_"$refSpikeIn"."
albam=$alprefix"Aligned.sortedByCoord.out.bam"
alsam=${albam::-3}"sam"
fsam=${alsam::-3}"filtered_vW1.sam"
fsamNSort=${fsam::-3}"Nsort.sam"
fbamRef=${fsam::-3}"Nsort.ref_allele.bam"
fbamAlt=${fsam::-3}"Nsort.alt_allele.bam"


# STAR alignment:
echo "Aligning..."

/home/amendelevich/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR \
     --readFilesIn $fastq1Path $fastq2Path \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix $alprefix \
     --runThreadN $nthreds \
     --outSAMtype BAM SortedByCoordinate \
     --varVCFfile $vcf \
     --outSAMattributes NH HI AS nM vA vG \
     --genomeDir $refDir \
     --outFilterMultimapNmax 1 \
     --waspOutputMode SAMtag \
     --quantMode GeneCounts

samtools view -h -o $alsam $albam

# Allelic reads filtering: 
head -1000 $alsam | grep "^@" > $fsam
grep "vW:i:1" $alsam >> $fsam

rm $alsam

# Allele separation:
echo "Separating alleles..."

samtools sort -n -o $fsamNSort $fsam

python3 /home/amendelevich/scripts/spike_in_proj/alleleseparation.py \
        --sam $fsamNSort --odir $alD \
        --paired 1 

rm $fsamNSort

for sam in $(ls $alD/*.sam); do
  bam=${sam::-3}bam
  samtools view -o $bam $sam
  rm $sam
done

# -------------
# Info filling:

star_log=$alprefix"Log.final.out"
assign_log=$alprefix"Aligned.sortedByCoord.out.filtered_vW1.Nsort.log_allelic_assignment.txt"
stats_out=$alprefix"alignment_stats.out"
chr_counts=$alprefix"chr_counts"

nFastq=$(grep "Number of input reads" $star_log | cut -f 2)
nAlign=$(grep "Uniquely mapped reads number" $star_log | cut -f 2)
pUMap=$(grep "Uniquely mapped reads %" $star_log | cut -f 2)
pMany=$(grep "% of reads mapped to too many loci" $star_log | cut -f 2)
pMiss=$(grep "% of reads unmapped: too many mismatches" $star_log | cut -f 2)
pShort=$(grep "% of reads unmapped: too short" $star_log | cut -f 2)
pOther=$(grep "% of reads unmapped: other" $star_log | cut -f 2)

assign_vect=($(head -9 $assign_log | tail -5 | cut -f1 -d " "))
nA1=${assign_vect[0]}
nA2=${assign_vect[1]}
nA1prob=${assign_vect[2]}
nA2prob=${assign_vect[3]}
nNoA=${assign_vect[4]}

echo -e $i"\t"$sampleID"\t"$refSample"\t"$refSpikeIn"\t"$nFastq"\t"$nAlign"\t"$pUMap"\t"$pMany"\t"$pMiss"\t"$pShort"\t"$pOther"\t"$albam"\t"$nA1"\t"$nA2"\t"$nA1prob"\t"$nA2prob"\t"$nNoA"\t"$fbamRef"\t"$fbamAlt 
echo -e $i"\t"$sampleID"\t"$refSample"\t"$refSpikeIn"\t"$nFastq"\t"$nAlign"\t"$pUMap"\t"$pMany"\t"$pMiss"\t"$pShort"\t"$pOther"\t"$albam"\t"$nA1"\t"$nA2"\t"$nA1prob"\t"$nA2prob"\t"$nNoA"\t"$fbamRef"\t"$fbamAlt > $stats_out

samtools view $albam   | cut -f3 | sort -V | uniq -c | sed -e 's/^ *//' -e 's/ /\t/' > $chr_counts".aligned.tsv"
samtools view $fbamRef | cut -f3 | sort -V | uniq -c | sed -e 's/^ *//' -e 's/ /\t/' > $chr_counts".ref_allele.tsv"
samtools view $fbamAlt | cut -f3 | sort -V | uniq -c | sed -e 's/^ *//' -e 's/ /\t/' > $chr_counts".alt_allele.tsv"


#!/bin/bash
#SBATCH -c 4
#SBATCH --mem-per-cpu=24G
#SBATCH -o /home/amendelevich/logs/bulkSpike/Sep22/fastq2albam_%A_%a.out
#SBATCH -e /home/amendelevich/logs/bulkSpike/Sep22/fastq2albam_%A_%a.err

module load gcc/5.3.0 
module load samtools/1.12 bcftools/1.7 htslib/1.7 python/3.6.10
module load R/3.6.1

# Table to fill:
infoTab=$1

indices=( $(cut -f1 $infoTab | tail -n +2) )

for i in ${indices[*]}; do

infoVect=(`grep -w "^"$i $infoTab`)

# 13 mandatory fields (0-12), so:
sampleID=${infoVect[1]}
refSample=${infoVect[4]}
refSpikeIn=${infoVect[5]}
refDir=${infoVect[9]}
fastq1Path=${infoVect[10]}
fastq2Path=${infoVect[11]}
d=${infoVect[12]}
suffChrSample=${infoVect[29]}
suffChrSpikein=${infoVect[30]}
 
vcf=$(ls $refDir/*vcf)
alD=$d/alignments/$sampleID/
alprefix=$alD/$sampleID".on_"$refSample"_"$refSpikeIn"."
albam=$alprefix"Aligned.sortedByCoord.out.bam"
alsam=${albam::-3}"sam"
fsam=${alsam::-3}"filtered_vW1.sam"
fsamNSort=${fsam::-3}"Nsort.sam"
fbamRef=${fsam::-3}"Nsort.ref_allele.bam"
fbamAlt=${fsam::-3}"Nsort.alt_allele.bam"

# -------------
# Info filling:

star_log=$alprefix"Log.final.out"
assign_log=$alprefix"Aligned.sortedByCoord.out.filtered_vW1.Nsort.log_allelic_assignment.txt"
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


chrcAling=$chr_counts".aligned.tsv"
chrcRef=$chr_counts".ref_allele.tsv"
chrcAlt=$chr_counts".alt_allele.tsv"

v29=$chr_counts".aligned_ref_alt.tsv"

script="df=Reduce(function(x,y) merge(x,y, by=\"V2\", all=T), list(data.frame(read.delim(\""$chrcAling"\", header=F)), data.frame(read.delim(\""$chrcRef"\", header=F)), data.frame(read.delim(\""$chrcAlt"\", header=F)))); write.table(df, file=\""$v29"\", sep=\"\t\", row.names=FALSE, col.names=F, quote=FALSE)"
R -e "$script"

function chromosomal_counter {
    precount=$(grep $1"$" $2 | cut -f1 | paste -sd+ | bc)
    [[ ! -z "$precount" ]] && echo $(($precount / 2)) || echo 0
}

v32=`chromosomal_counter $suffChrSample $chrcAling`
v33=`chromosomal_counter $suffChrSpikein $chrcAling`
v34=`chromosomal_counter $suffChrSample $chrcRef`
v35=`chromosomal_counter $suffChrSample $chrcAlt`
v36=`chromosomal_counter $suffChrSpikein $chrcRef`
v37=`chromosomal_counter $suffChrSpikein $chrcAlt`

#v32=$(($(grep $suffChrSample"$" $chrcAling | cut -f1 | paste -sd+ | bc) / 2))
#v33=$(($(grep $suffChrSpikein"$" $chrcAling | cut -f1 | paste -sd+ | bc) / 2))
#v34=$(($(grep $suffChrSample"$" $chrcRef | cut -f1 | paste -sd+ | bc) / 2))
#v35=$(($(grep $suffChrSample"$" $chrcAlt | cut -f1 | paste -sd+ | bc) / 2))
#v36=$(($(grep $suffChrSpikein"$" $chrcRef | cut -f1 | paste -sd+ | bc) / 2))
#v37=$(($(grep $suffChrSpikein"$" $chrcAlt | cut -f1 | paste -sd+ | bc) / 2))

# columns 14-28, 29 and 32-37 

awk -v id=$sampleID -v d=$d \
    -v v14=$nFastq -v v15=$nAlign -v v16=$pUMap \
    -v v17=$pMany   -v v18=$pMiss  -v v19=$pShort -v v20=$pOther \
    -v v21=$albam \
    -v v22=$nA1    -v v23=$nA2    -v v24=$nA1prob -v v25=$nA2prob -v v26=$nNoA \
    -v v27=$fbamRef -v v28=$fbamAlt \
    -v v29=$v29 \
    -v v32=$v32 -v v33=$v33 -v v34=$v34 -v v35=$v35 -v v36=$v36 -v v37=$v37 \
    -v OFS='\t' '{
        if ($2==id && $13==d) {
            $14=v14; $15=v15; $16=v16; $17=v17; $18=v18; $19=v19; $20=v20; $21=v21;
            $22=v22; $23=v23; $24=v24; $25=v25; $26=v26; $27=v27; $28=v28;
            $29=v29; $32=v32; $33=v33; $34=v34; $35=v35; $36=v36; $37=v37
        }; print $0}' $infoTab > $d/tmp$id; \
    mv $d/tmp$id $infoTab 

done

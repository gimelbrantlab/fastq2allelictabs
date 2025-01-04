#!/bin/bash

# EXAMPLE:  
#   bash $scripts_d/fastq2allelicbams_stats.sh $TAB 
#   bash ~/my_dir/fastq2allelicbams_stats.sh /path/to/table_info.tsv 

# Tables (input and output):
infoTab=$1
infoTab_filled=${infoTab::-4}_info.tsv

# Indices from the 1st column:
indices=( $(cut -f1 $infoTab | tail -n +2) )

# New column names (only 1-15 are mandatory be present in the initial table):
colnames=('##' sampleID platform readLenPE refSample refSpikeIn propSpikeIn chimericName refDirA1 refDirA2 fastq1Path fastq2Path outPath suffChrSample suffChrSpikein bamAlignAR1 bamAlignAR2 bamA1 bamA2 tabNChrom nFastq nAlignAR1 nAlignAR1Sample nAlignAR1Spikein pUMapAR1 pManyAR1 pMissAR1 pShortAR1 pOtherAR1 nAlignAR2 nAlignAR2Sample nAlignAR2Spikein pUMapAR2 pManyAR2 pMissAR2 pShortAR2 pOtherAR2 nA1 nA2 nNoA nA1Sample nA2Sample nA1Spikein nA2Spikein)

# ------ Run stats collection row by row ---------------------------------------------------------------------
for i in ${indices[*]}; do # ---------------------------------------------------------------------------------

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
suffChrSample=${infoVect[13]}
suffChrSpikein=${infoVect[14]}
 
vcf=$(ls $refDir1/*vcf)
alD=$d/alignments/$sampleID/

# -- computing stats for alignments on A1 and A2 --
# --
# -- files we need --

alprefix=$alD/$sampleID".on_"$chimericName
sampleprefix=$alD/$sampleID"_"$chimericName

bamAlignAR1=$alprefix".on_Allele1.Aligned.sortedByCoord.out.bam"
bamAlignAR2=$alprefix".on_Allele2.Aligned.sortedByCoord.out.bam"

bamA1=$sampleprefix".Allele1.bam"
bamA2=$sampleprefix".Allele2.bam"

star_log_AR1=$alprefix".on_Allele1.Log.final.out"
star_log_AR2=$alprefix".on_Allele2.Log.final.out"

assign_log=$sampleprefix".log_allelic_assignment.txt"

chr_counts=$sampleprefix".chr_counts"
chrcAlign1=$chr_counts".aligned_on_Allele1.tsv"
chrcAlign2=$chr_counts".aligned_on_Allele2.tsv"
chrcRef=$chr_counts".Allele1.tsv"
chrcAlt=$chr_counts".Allele2.tsv"

# -- taking out all stats --
# ---- STAR log -- 

nFastq=$(grep "Number of input reads" $star_log_AR1 | cut -f 2)

nAlignAR1=$(grep "Uniquely mapped reads number" $star_log_AR1 | cut -f 2)
pUMapAR1=$(grep "Uniquely mapped reads %" $star_log_AR1 | cut -f 2)
pManyAR1=$(grep "% of reads mapped to too many loci" $star_log_AR1 | cut -f 2)
pMissAR1=$(grep "% of reads unmapped: too many mismatches" $star_log_AR1 | cut -f 2)
pShortAR1=$(grep "% of reads unmapped: too short" $star_log_AR1 | cut -f 2)
pOtherAR1=$(grep "% of reads unmapped: other" $star_log_AR1 | cut -f 2)

nAlignAR2=$(grep "Uniquely mapped reads number" $star_log_AR2 | cut -f 2)
pUMapAR2=$(grep "Uniquely mapped reads %" $star_log_AR2 | cut -f 2)
pManyAR2=$(grep "% of reads mapped to too many loci" $star_log_AR2 | cut -f 2)
pMissAR2=$(grep "% of reads unmapped: too many mismatches" $star_log_AR2 | cut -f 2)
pShortAR2=$(grep "% of reads unmapped: too short" $star_log_AR2 | cut -f 2)
pOtherAR2=$(grep "% of reads unmapped: other" $star_log_AR2 | cut -f 2)


# ---- merge log --
asums=(`tail -n +5 $assign_log | head -n 15 | cut -f1 -d " "`)
nA1=$(( ${asums[0]} + ${asums[2]} ))
nA2=$(( ${asums[1]} + ${asums[3]} ))
nNoA=${asums[4]}

# ---- sample and spikein log --

# ------ make the a1 and a2 merged table of chromosomes counts
tabNChrom=$chr_counts".aligned_ref_alt.tsv"

script="df=Reduce(function(x,y) merge(x,y, by=\"V2\", all=T), list(data.frame(read.delim(\""$chrcAlign1"\", header=F)), data.frame(read.delim(\""$chrcAlign2"\", header=F)), data.frame(read.delim(\""$chrcRef"\", header=F)), data.frame(read.delim(\""$chrcAlt"\", header=F)))); write.table(df, file=\""$tabNChrom"\", sep=\"\t\", row.names=FALSE, col.names=F, quote=FALSE)"
R -e "$script"

# ------ count sample and spike-in coverage separately
function chromosomal_counter {
    if [[ "$1" == "None" ]]; then
        echo 0 ; return  
        # no spike-in chromosomes to align on them
    elif [[ "$1" == "-" ]]; then
        precount=$(cut -f1 $2 | awk '{sum += $1} END {print sum}')  
        # only sample chrs in reference -> sum up everything
    else
        precount=$(grep $1"$" $2 | cut -f1 | awk '{sum += $1} END {print sum}')
    fi
    [[ ! -z "$precount" ]] && echo $(($precount / 2)) || echo 0
                                              # /2 because it's PE, and we're counting fragments
}

nAlignAR1Sample=`chromosomal_counter $suffChrSample $chrcAlign1`
nAlignAR1Spikein=`chromosomal_counter $suffChrSpikein $chrcAlign1`
nAlignAR2Sample=`chromosomal_counter $suffChrSample $chrcAlign2`
nAlignAR2Spikein=`chromosomal_counter $suffChrSpikein $chrcAlign2`
nA1Sample=`chromosomal_counter $suffChrSample $chrcRef`
nA2Sample=`chromosomal_counter $suffChrSample $chrcAlt`
nA1Spikein=`chromosomal_counter $suffChrSpikein $chrcRef`
nA2Spikein=`chromosomal_counter $suffChrSpikein $chrcAlt`

# -- writing down -----------
# -- filling columns 16-44 --

awk -v id=$sampleID -v d=$d \
    -v v16=$bamAlignAR1 -v v17=$bamAlignAR2 \
    -v v18=$bamA1 -v v19=$bamA2 \
    -v v20=$tabNChrom \
    -v v21=$nFastq \
    -v v22=$nAlignAR1 -v v23=$nAlignAR1Sample -v v24=$nAlignAR1Spikein \
    -v v25=$pUMapAR1  -v v26=$pManyAR1        -v v27=$pMissAR1         \
    -v v28=$pShortAR1 -v v29=$pOtherAR1                                \
    -v v30=$nAlignAR2 -v v31=$nAlignAR2Sample -v v32=$nAlignAR2Spikein \
    -v v33=$pUMapAR2  -v v34=$pManyAR2        -v v35=$pMissAR2         \
    -v v36=$pShortAR2 -v v37=$pOtherAR2                                \
    -v v38=$nA1       -v v39=$nA2             -v v40=$nNoA             \
    -v v41=$nA1Sample -v v42=$nA2Sample -v v43=$nA1Spikein -v v44=$nA2Spikein \
    -v OFS='\t' '{
        if ($2==id && $13==d) {
            $16=v16; $17=v17; $18=v18; $19=v19; $20=v20; $21=v21;
            $22=v22; $23=v23; $24=v24; $25=v25; $26=v26; $27=v27; $28=v28;
            $29=v29; $30=v30; $31=v31; $32=v32; $33=v33; $34=v34; $35=v35; 
            $36=v36; $37=v37; $38=v38; $39=v39; $40=v40; $41=v41; $42=v42;
            $43=v43; $44=v44 
        }; print $0}' $infoTab_filled > $d/tmp$id; \
    
    (echo "${colnames[*]}" | tr ' ' '\t'; tail -n +2 $d/tmp$id) > $infoTab_filled

    rm $d/tmp$id 

done # ---------------------------------------------------------------------------------

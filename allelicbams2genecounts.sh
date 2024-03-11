#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu=24G
#SBATCH -o /home/amendelevich/logs/palladium/13apr23/genecount_%A_%a.out
#SBATCH -e /home/amendelevich/logs/palladium/13apr23/genecount_%A_%a.err

module load gcc samtools

# Info Table:
infoTab=$1

infoName=$(basename $infoTab)
outName=${infoName::-9}

info_uniq=${infoTab::-4}".uniq4featurecounts.tsv"
cut -f5-6,13 $infoTab | tail -n +2 | sort -u > $info_uniq

while read line; do
    xline=(${line[0]}) 
    org_sample=${xline[0]}
    org_spikein=${xline[1]}
    ddir=${xline[2]}
    odir=$ddir/featurecounts
    mkdir $odir
    otxt=$odir/$outName"."$org_sample"_"$org_spikein".counts.tsv"
    bams=($( awk -v os1=$org_sample -v os2=$org_spikein -v d=$ddir -v OFS='\t' '{if ($5==os1 && $6==os2 && $13==d) {print $18, $19}}' $infoTab ))
    echo ${bams[*]} "\n"
    refdirs=($( awk -v os1=$org_sample -v os2=$org_spikein -v d=$ddir -v OFS='\t' '{if ($5==os1 && $6==os2 && $13==d) {print $10}}' $infoTab ))

    echo /home/amendelevich/tools/subread-2.0.2-Linux-x86_64/bin/featureCounts -p -a $(ls ${refdirs[0]}/*gtf) -o $otxt ${bams[*]} --countReadPairs -B -C "\n"

    /home/amendelevich/tools/subread-2.0.2-Linux-x86_64/bin/featureCounts -p -a $(ls ${refdirs[0]}/*gtf) -o $otxt ${bams[*]} --countReadPairs -B -C

    ### same for original alignments (on allele1, on allele2):
    
    otxt=$odir/$outName"."$org_sample"_"$org_spikein".total_counts.ref_A1.tsv"
    bams=($( awk -v os1=$org_sample -v os2=$org_spikein -v d=$ddir -v OFS='\t' '{if ($5==os1 && $6==os2 && $13==d) {print $16}}' $infoTab ))
    echo ${bams[*]} "\n"
    echo /home/amendelevich/tools/subread-2.0.2-Linux-x86_64/bin/featureCounts -p -a $(ls ${refdirs[0]}/*gtf) -o $otxt ${bams[*]} --countReadPairs -B -C "\n"
    /home/amendelevich/tools/subread-2.0.2-Linux-x86_64/bin/featureCounts -p -a $(ls ${refdirs[0]}/*gtf) -o $otxt ${bams[*]} --countReadPairs -B -C

    otxt=$odir/$outName"."$org_sample"_"$org_spikein".total_counts.ref_A2.tsv"
    bams=($( awk -v os1=$org_sample -v os2=$org_spikein -v d=$ddir -v OFS='\t' '{if ($5==os1 && $6==os2 && $13==d) {print $17}}' $infoTab ))
    echo ${bams[*]} "\n"
    echo /home/amendelevich/tools/subread-2.0.2-Linux-x86_64/bin/featureCounts -p -a $(ls ${refdirs[0]}/*gtf) -o $otxt ${bams[*]} --countReadPairs -B -C "\n"
    /home/amendelevich/tools/subread-2.0.2-Linux-x86_64/bin/featureCounts -p -a $(ls ${refdirs[0]}/*gtf) -o $otxt ${bams[*]} --countReadPairs -B -C

    ### same for undetermined alignments:
    ###################################

done < $info_uniq








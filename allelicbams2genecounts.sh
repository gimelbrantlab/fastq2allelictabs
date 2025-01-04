#!/bin/bash

# Example:
#   bash $scripts_d/allelicbams2genecounts.sh ${TAB_with_stats}  
# Output: featurecounts tables for all unique pairs of species. 

# Info Table with stats added on the previous point (so its name is xxx_info.tsv):
infoTab=$1

infoName=$(basename $infoTab)
outName=${infoName::-9}

# Form groups from unique pairs of species (sample and spike-in):

info_uniq=${infoTab::-4}".uniq4featurecounts.tsv"
cut -f5-6,13 $infoTab | tail -n +2 | sort -u > $info_uniq

# Run featurecounts for each group:

while read line; do
    xline=(${line[0]}) 
    org_sample=${xline[0]}
    org_spikein=${xline[1]}
    ddir=${xline[2]}
    odir=$ddir/featurecounts
    mkdir $odir
    otxt=$odir/$outName"."$org_sample"_"$org_spikein".counts.tsv"
    
    # select bams and fastas for a group:
    bams=($( awk -v os1=$org_sample -v os2=$org_spikein -v d=$ddir -v OFS='\t' '{if ($5==os1 && $6==os2 && $13==d) {print $18, $19}}' $infoTab ))
    echo ${bams[*]} "\n"
    refdirs=($( awk -v os1=$org_sample -v os2=$org_spikein -v d=$ddir -v OFS='\t' '{if ($5==os1 && $6==os2 && $13==d) {print $10}}' $infoTab ))

    # run featurecouns:
    featureCounts -p -a $(ls ${refdirs[0]}/*gtf) -o $otxt ${bams[*]} --countReadPairs -B -C 

done < $info_uniq








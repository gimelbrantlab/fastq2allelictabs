#########################################################
#                                                       #
#                 ASE Spike-In Project                  #
#                                                       #
#      Authors: Mendelevich Asia, Pakharev Aleksei	#
#                                                       #
#           Earlier versions: ASEReadcounter*           #
#               with Svetlana Vinogradova               #
#     Idea is originaly taken from Stephane Castel      #
#                                                       #
#########################################################
##
## NOTE: 
## 1. SAM-files should be already sorted by read names
## 2. Read names must be non-empty (to preserve the best fitting reference) 
## 3. AS field must be present in SAM files (to filter by alignment quality) - typical for STAR
## 4. vA field must be present in SAM files (to filter by heterozygous SNPs) - the respective flag should be used in STAR
## 
## USAGE:  python3 alleleseparation.py --samA1 [path to alignment on A1] --samA2 [path to alignment on A2] --obase [meaningful part of the name, used in outputs] --odir [output dir] --allele1 [name of the first allele] --allele2 [same for allele 2] --paired [0|1] 
##

#import cProfile
import sys
import pandas
import argparse
import random
import subprocess
import time
import os
start_time = time.time()

AS_prefix = "AS:i:" 
# quality, accounting for mismatches
vA_prefix = "vA:B:c,"
# from STAR manual: 1 or 2 match one of the genotype alleles, 3 - no match to genotype.
vG_prefix = "vG:B:i,"
# SNP positions (without chr, but it's well-determined by pos only)

def get_name_score_variants(line):
    name_pos = line.find("\t")
    score_start_pos = line.find(AS_prefix, name_pos) + len(AS_prefix)
    score_fin_pos = line.find("\t", score_start_pos)
    vA_start_pos = line.find(vA_prefix, name_pos) + len(vA_prefix)
    vA_end_pos   = line.find('\t', vA_start_pos)
    vG_start_pos = line.find(vG_prefix, name_pos) + len(vG_prefix)
    vG_end_pos   = line.find('\t', vG_start_pos)
    return line[:name_pos], int(line[score_start_pos:score_fin_pos]), line[vA_start_pos:vA_end_pos], line[vG_start_pos:vG_end_pos]

def parse_variants(variants, positions):
    # rule1: % of minor allele (ref/alt) is <= 15 %
    # rule2: % of major allele (ref/alt) is >= 80 %
    vectvars = variants.strip().split(",")
    vectpos  = positions.strip().split(",")
    ### is reads intersect and len(fragm)<2*len(read),
    ### then it might have same snp counted twice,
    ### we want to take unique encountersm and remove reads
    ### having mismatching calls at snps:

    dictvars = dict()
    clean = True
    for i, pos in enumerate(vectpos):
        if (pos in dictvars) and (dictvars[pos] != vectvars[i]):
            clean = False
        else:
            dictvars[pos] = vectvars[i]

    l = len(dictvars)
    if not clean:
        return ("disc", 0,  {"c_sum":l, "c_this":0})

    c1 , c2 , c3 = [list(dictvars.values()).count(x) for x in ['1','2','3']]
    if c1 == l:
        return ("only" , 1, {"c_sum":l, "c_this":l})
    elif c2 == l:
        return ("only" , 2, {"c_sum":l, "c_this":l})
    elif c1/l >= 0.8 and c2/l <= 0.15:
        return ("best" , 1, {"c_sum":l, "c_this":c1})
    elif c2/l >= 0.8 and c1/l <= 0.15: 
        return ("best" , 2, {"c_sum":l, "c_this":c2})
    else:
        return ("none" , 0, {"c_sum":l, "c_this":0})

# Dealing with samtools ordering by name:

def isdigit(c):
    return c.isdigit()

def uporotiycmp(firstrb, secstrb):
    firstr = firstrb + '\x00' 
    secstr = secstrb + '\x00'
    firind , secind = 0 , 0
    firlen = len(firstr)
    seclen = len(secstr)
    while (firind < firlen) and (secind < seclen):
        if (isdigit(firstr[firind]) and isdigit(secstr[secind])):
            while (firstr[firind] == '0'):
                firind += 1
            while (secstr[secind] == '0'):
                secind += 1
            while isdigit(firstr[firind]) and isdigit(secstr[secind]) and (firstr[firind] == secstr[secind]):
                firind += 1
                secind += 1
            if (isdigit(firstr[firind]) and isdigit(secstr[secind])):
                i = 0
                while (isdigit(firstr[firind+i]) and isdigit(secstr[secind+i])):
                    i += 1
                if isdigit(firstr[firind+i]): return 1
                elif isdigit(secstr[secind+i]): return -1
                else: return ord(firstr[firind]) - ord(secstr[secind])
            elif (isdigit(firstr[firind])): return 1
            elif (isdigit(secstr[secind])): return -1
            elif (firind < secind): return 1
            elif (firind > secind): return -1
        else:
            if (firstr[firind] != secstr[secind]):
               return ord(firstr[firind]) - ord(secstr[secind])
            firind += 1
            secind += 1
    if firind < firlen: return 1
    elif secind < seclen: return -1
    else: return 0

def asc_order(name1, name2):
    return (uporotiycmp(name1, name2) <= 0)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samA1", required=True, help="Path to reads aligned to Allele1 pseudogenome")
    parser.add_argument("--samA2", required=True, help="Path to reads aligned to Allele2 pseudogenome")
    parser.add_argument("--odir", required=True, help="Path to output directory")
    parser.add_argument("--obase", required=True, help="Output base name (prefix for outputs)")
    parser.add_argument("--allele1", default="ref_allele", help="Name of the first allele (default: ref_allele)")
    parser.add_argument("--allele2", default="alt_allele", help="Name of the second allele (default: alt_allele)")
    parser.add_argument("--paired", default=1, help="Flag: If reads are paired-end (yes : 1, no : 0)")
    args = parser.parse_args()
    
    paired = int(args.paired)

    # Output files will be (confusing namings, yep):

    output_name_base = args.obase
    sam1a = os.path.join(args.odir, ".".join([output_name_base, args.allele1, "sam"]))
    sam2a = os.path.join(args.odir, ".".join([output_name_base, args.allele2, "sam"]))
    samNa = os.path.join(args.odir, ".".join([output_name_base, "ambiguous", "sam"]))

    # Open output sams; Get header:

    header = subprocess.check_output("samtools view -SH --no-PG " + args.samA1, shell=True, universal_newlines=True)
    rg_a2  = subprocess.check_output("samtools view -SH --no-PG " + args.samA2 + '| grep "^@RG"', shell=True, universal_newlines=True)

    out_stream_1a = open(sam1a, "w")
    out_stream_2a = open(sam2a, "w")
    out_stream_Na = open(samNa, "w")
    out_stream_1a.write(header)
    out_stream_2a.write(header)
    out_stream_Na.write(header)
    out_stream_1a.write(rg_a2)
    out_stream_2a.write(rg_a2)
    out_stream_Na.write(rg_a2)

    # Open input_sam:        

    a1_counts = {"clearSNPs":0, "singleAllele":0, "singleAllele_clearSNPs":0, "rest":0}
    a2_counts = {"clearSNPs":0, "singleAllele":0, "singleAllele_clearSNPs":0, "rest":0}
    nondetermined_count = {"disc_mates":0 , "disc_snp":0 , "u_disc_snp_a1":0 , "u_disc_snp_a2":0 , "diff_position":0 , "disc_snp_a1_q_a2":0 , "disc_snp_a2_q_a1":0}
    bad_reads = set()

    err_2_disc_classes = ["disc_snp", "disc_snp_a1_q_a2", "disc_snp_a2_q_a1"]

    source_A1 = open(args.samA1, 'r')
    source_A2 = open(args.samA2, 'r')

    # Skip header in each file:
    
    A1_skip = int(subprocess.check_output("samtools view -SH --no-PG "+args.samA1+" | wc -l", shell=True, universal_newlines=True).strip())
    A2_skip = int(subprocess.check_output("samtools view -SH --no-PG "+args.samA2+" | wc -l", shell=True, universal_newlines=True).strip())

    for i in range(A1_skip):
        source_A1.readline()
    for i in range(A2_skip):
        source_A2.readline()

    # Read process:

    def output_read(fhandler, a_read):
        if paired:
            fhandler.write(a_read[0])
            fhandler.write(a_read[1])
        else:
            fhandler.write(a_read)

    def blocks_generator(fhandler):
        beg_line = fhandler.readline()
        output = [beg_line]
        name, score, vA, vG = get_name_score_variants(beg_line)
        our_prefix = name+'\t'
        for line in fhandler:
            if line.startswith(our_prefix):
                output.append(line)
            else:
                yield output, name, score, vA, vG
                beg_line = line
                output = [beg_line]
                name, score, vA, vG = get_name_score_variants(beg_line)
                our_prefix = name+'\t'
        yield output, name, score, vA, vG

    def correct_blocks_generator(fhandler):
        if paired:
            for block, name, score, vA, vG in blocks_generator(fhandler):
                if (len(block) == 2): yield block, name, score, vA, vG
                else: bad_reads.add(name)
        else:
            for block, name, score, vA, vG in blocks_generator(fhandler):
                if (len(block) == 1): yield block[0], name, score, vA, vG
                else: bad_reads.add(name)

    # Create generator objects:
    
    A2gen = correct_blocks_generator(source_A2)
    A1gen = correct_blocks_generator(source_A1)
    
  
    # Separate till any EOF:
    try:
        a1_read, a1_read_name, a1_score, a1_vA, a1_vG = A1gen.__next__()
        a2_read, a2_read_name, a2_score, a2_vA, a2_vG = A2gen.__next__()
        while 1:
            if (a1_read_name == a2_read_name):
                # i.e aligned to both alleles
                # are they aligned to the same place (SNPs)?:
                if (a1_vG != a2_vG or a1_vA != a2_vA):
                    nondetermined_count["diff_position"] += 1
                    output_read(out_stream_Na, a1_read)
                    output_read(out_stream_Na, a2_read)
                else:
                    vA1 = parse_variants(a1_vA, a1_vG)
                    vA2 = parse_variants(a2_vA, a2_vG)
                    # can this read be assigned to some allele looking at SNPs?
                    # is the best candidate when looking at SNPs has better quality too?
                    if (vA1[1] == 1 and (a1_score >= a2_score or a1_score <= a2_score and vA1[2]["c_this"] >= 3)):
                        if (vA1[0] == "only"):
                            a1_counts["clearSNPs"] += 1
                        else:
                            a1_counts["rest"] += 1
                        output_read(out_stream_1a, a1_read)
                    elif (vA2[1] == 2 and (a1_score <= a2_score or a1_score >= a2_score and vA2[2]["c_this"] >= 3)):
                        if (vA2[0] == "only"):
                            a2_counts["clearSNPs"] += 1
                        else:
                            a2_counts["rest"] += 1
                        output_read(out_stream_2a, a2_read)
                    elif (vA1[0] == "disc"):
                        nondetermined_count["disc_mates"] += 1
                        output_read(out_stream_Na, a1_read)
                        output_read(out_stream_Na, a2_read)
                    else:
                        nondetermined_count[err_2_disc_classes[vA1[1]]] += 1
                        output_read(out_stream_Na, a1_read)
                        output_read(out_stream_Na, a2_read)
                a1_read, a1_read_name, a1_score, a1_vA, a1_vG = A1gen.__next__() 
                a2_read, a2_read_name, a2_score, a2_vA, a2_vG = A2gen.__next__()            
            elif asc_order(a1_read_name, a2_read_name):
                # i.e should look at a1_read now (aligned to A1 only)
                vA = parse_variants(a1_vA, a1_vG)
                if (vA[1] == 1):
                    if(vA[0] == "only"):
                        a1_counts["singleAllele_clearSNPs"] += 1
                    else:
                        a1_counts["singleAllele"] += 1
                    sout = out_stream_1a
                elif (vA[0] == "disc"):
                    nondetermined_count["disc_mates"] += 1
                    sout = out_stream_Na
                else:
                    nondetermined_count["u_disc_snp_a1"] += 1
                    sout = out_stream_Na
                output_read(sout, a1_read)
                a1_read, a1_read_name , a1_score , a1_vA , a1_vG = A1gen.__next__()
            else:
                # i.e should look at a2_read now (aligned to A2 only)
                vA = parse_variants(a2_vA, a2_vG)
                if (vA[1] == 2):
                    if(vA[0] == "only"):
                        a2_counts["singleAllele_clearSNPs"] += 1
                    else:
                        a2_counts["singleAllele"] += 1
                    sout = out_stream_2a
                elif (vA[0] == "disc"):
                    nondetermined_count["disc_mates"] += 1
                    sout = out_stream_Na
                else:
                    nondetermined_count["u_disc_snp_a2"] += 1
                    sout = out_stream_Na
                output_read(sout, a2_read)
                a2_read, a2_read_name , a2_score , a2_vA , a2_vG = A2gen.__next__() 
    except StopIteration:
        pass

    # ... and write the remain part of reads, checking vA consistence with allele:
    for a1_read, a1_read_name , a1_score , a1_vA , a1_vG in A1gen:
        vA = parse_variants(a1_vA, a1_vG)
        if (vA[1] == 1):
            if(vA[0] == "only"):
                a1_counts["singleAllele_clearSNPs"] += 1
            else:
                a1_counts["singleAllele"] += 1
            sout = out_stream_1a
        elif (vA[0] == "disc"):
            nondetermined_count["disc_mates"] += 1
            sout = out_stream_Na
        else:
            nondetermined_count["u_disc_snp_a1"] += 1
            sout = out_stream_Na
        output_read(sout, a1_read)
    for a2_read, a2_read_name , a2_score , a2_vA , a2_vG in A2gen:
        vA = parse_variants(a2_vA, a2_vG)
        if (vA[1] == 2):
            if(vA[0] == "only"):
                a2_counts["singleAllele_clearSNPs"] += 1
            else:
                a2_counts["singleAllele"] += 1
            sout = out_stream_2a
        elif (vA[0] == "disc"):
            nondetermined_count["disc_mates"] += 1
            sout = out_stream_Na
        else:
            nondetermined_count["u_disc_snp_a2"] += 1
            sout = out_stream_Na
        output_read(sout, a2_read)


    # Closing:

    out_stream_Na.close(); out_stream_1a.close(); out_stream_2a.close()
    source_A1.close(); source_A2.close()

    logfile = os.path.join(args.odir, ".".join([output_name_base,"log_allelic_assignment", "txt"]))
    logtime = time.time() - start_time

    outmessage = []
    outmessage.append("----------------------------------------------------------------------------------------")
    outmessage.append("ALLELIC ASSIGNMENT SUMMARY for \n%s"%(output_name_base))
    outmessage.append("----------------------------------------------------------------------------------------")
    outmessage.append("%d reads only in alignment on Allele 1"%(a1_counts["singleAllele_clearSNPs"]+a1_counts["singleAllele"]))
    outmessage.append("    %d - all SNPs"%(a1_counts["singleAllele_clearSNPs"]))
    outmessage.append("    %d - dominating number of SNPs"%(a1_counts["singleAllele"]))
    outmessage.append("%d reads only in alignment on Allele 2"%(a2_counts["singleAllele_clearSNPs"]+a2_counts["singleAllele"]))
    outmessage.append("    %d - all SNPs"%(a2_counts["singleAllele_clearSNPs"]))
    outmessage.append("    %d - dominating number of SNPs"%(a2_counts["singleAllele"]))
    outmessage.append("%d reads better aligned to Allele 1"%(a1_counts["clearSNPs"]+a1_counts["rest"]))
    outmessage.append("    %d - all SNPs"%(a1_counts["clearSNPs"]))
    outmessage.append("    %d - dominating number of SNPs"%(a1_counts["rest"]))
    outmessage.append("%d reads better aligned to Allele 2"%(a2_counts["clearSNPs"]+a2_counts["rest"]))
    outmessage.append("    %d - all SNPs"%(a2_counts["clearSNPs"]))
    outmessage.append("    %d - dominating number of SNPs"%(a2_counts["rest"]))
    outmessage.append("%d reads were not determined"%(sum(nondetermined_count.values())))
    outmessage.append("    %d - unmatching alignment positioning between alleles"%(nondetermined_count["diff_position"]))
    outmessage.append("    %d - SNPs vote for A1, quality votes for A2"%(nondetermined_count["disc_snp_a1_q_a2"]))
    outmessage.append("    %d - SNPs vote for A2, quality votes for A1"%(nondetermined_count["disc_snp_a2_q_a1"]))
    outmessage.append("    %d - was aligned only to A1, but SNPs come from A2"%(nondetermined_count["u_disc_snp_a1"]))
    outmessage.append("    %d - was aligned only to A2, but SNPs come from A1"%(nondetermined_count["u_disc_snp_a2"]))
    outmessage.append("    %d - discordant SNP calls within fragment"%(nondetermined_count["disc_snp"]))
    outmessage.append("    %d - discordant SNP calls on mate reads"%(nondetermined_count["disc_mates"]))
    outmessage.append("----- %s seconds -----" %(logtime))
    outmessage.append("%d BAD read names: "%(len(bad_reads)) + " , ".join(sorted(list(bad_reads))))
    outmessage.append("----------------------------------------------------------------------------------------\n")
    print("\n".join(outmessage))

    with open (logfile, "w") as f:
        f.write("\n".join(outmessage))


if __name__ == "__main__":
    #cProfile.run("main()")
    main()

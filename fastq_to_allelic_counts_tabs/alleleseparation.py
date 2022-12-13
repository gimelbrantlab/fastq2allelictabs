#########################################################
#                                                       #
#                 ASE Spike-In Project                  #
#                                                       #
#               Author: Mendelevich Asia		#
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
## USAGE:  python3 alleleseparation.py --samA1 [path to alignment on A1] --samA2 [path to alignment on A2] --odir [output dir] --allele1 [name of the first allele] --allele2 [same for allele 2] --paired [0|1] 
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

vA_prefix = "vA:B:c,"
# from STAR manual: 1 or 2 match one of the genotype alleles, 3 - no match to genotype.

def get_name_and_variants(line):
    vA_start_pos = line.find(vA_prefix) + len(vA_prefix)
    vA_end_pos   = line.find('\t', vA_start_pos)
    return line[:line.find('\t')], line[vA_start_pos:vA_end_pos]

def parse_variants(variants):
    # rule1: % of minor allele (ref/alt) is <= 15 %
    # rule2: % of major allele (ref/alt) is >= 80 %
    vectvars = variants.strip().split(",")
    l = len(vectvars)
    c1 , c2 , c3 = [vectvars.count(x) for x in ['1','2','3']]
    if c1 == l:
        return ("only" , 1)
    elif c2 == l:
        return ("only" , 2)
    elif c1/l >= 0.8 and c2/l <= 0.15:
        return ("best" , 1)
    elif c2/l >= 0.8 and c1/l <= 0.15: 
        return ("best" , 2)
    else:
        return ("none" , 0)



def main():
    # Parse arguments:
    # --sam [path to file] --odir [output dir] --allele1 [name of the first allele] --allele2 [same for allele 2] --paired [0|1] 

    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", required=True, help="Path to reads aligned with allelic STAR version, to allele1 genome and with provided heterozygous SNPs (vcf file with ref=allele1, alt=allele2)")
    parser.add_argument("--odir", required=True, help="Path to output directory")
    parser.add_argument("--allele1", default="ref_allele", help="Name of the first allele (default: ref_allele)")
    parser.add_argument("--allele2", default="alt_allele", help="Name of the second allele (default: alt_allele)")
    parser.add_argument("--paired", default=0, help="Flag: If reads are paired-end (yes : 1, no : 0)")
    args = parser.parse_args()
    
    paired = int(args.paired)

    output_name_base = os.path.basename(args.sam)[:-4]
    sam1a = os.path.join(args.odir, ".".join([output_name_base, args.allele1, "sam"]))
    sam2a = os.path.join(args.odir, ".".join([output_name_base, args.allele2, "sam"]))
    samNa = os.path.join(args.odir, ".".join([output_name_base, "ambiguous", "sam"]))

    # Open output sams; Get header:

    header = subprocess.check_output("samtools view -SH "+args.sam, shell=True, universal_newlines=True)
    out_stream_1a = open(sam1a, "w")
    out_stream_2a = open(sam2a, "w")
    out_stream_Na = open(samNa, "w")
    out_stream_1a.write(header)
    out_stream_2a.write(header)
    out_stream_Na.write(header)

    # Open input_sam:        

    a1_only , a1_count , a2_only , a2_count , nondetermined_count = 0 , 0 , 0 , 0 , 0
    bad_reads = set()

    source_s = open(args.sam, 'r')

    # Skip header in the input sam file:

    s_skip = int(subprocess.check_output("samtools view -SH "+args.sam+" | wc -l", shell=True, universal_newlines=True).strip())

    for i in range(s_skip-1):
        source_s.readline()
    
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
        name, vA = get_name_and_variants(beg_line)
        our_prefix = name+'\t'
        for line in fhandler:
            if line.startswith(our_prefix):
                output.append(line)
            else:
                yield output, name, vA
                beg_line = line
                output = [beg_line]
                name, vA = get_name_and_variants(beg_line)
                our_prefix = name+'\t'
        yield output, name, vA

    def correct_blocks_generator(fhandler):
        if paired:
            for block, name, vA in blocks_generator(fhandler):
                if (len(block) == 2): yield block, name, vA
                else: bad_reads.add(name)
        else:
            for block, name, vA in blocks_generator(fhandler):
                if (len(block) == 1): yield block[0], name, vA
                else: bad_reads.add(name)

    # Create generator objects:
    
    sgen = correct_blocks_generator(source_s)
      
    # Separate till EOF:
    try:
        s_read, s_readname, s_vA = sgen.__next__()
        while 1:
            vA = parse_variants(s_vA)
            if vA[1] == 1:
                if vA[0] == "only":
                    a1_only += 1
                elif vA[0] == "best":
                    a1_count += 1
                sout = out_stream_1a
            elif vA[1] == 2:
                if vA[0] == "only":
                    a2_only += 1
                elif vA[0] == "best":
                    a2_count += 1
                sout = out_stream_2a
            else:
                nondetermined_count += 1
                sout = out_stream_Na
                
            output_read(sout, s_read)
            s_read, s_readname, s_vA = sgen.__next__()
    except StopIteration:
        pass

    # Closing:

    out_stream_Na.close(); out_stream_1a.close(); out_stream_2a.close(); source_s.close()

    logfile = os.path.join(args.odir, ".".join([output_name_base,"log_allelic_assignment", "txt"]))
    logtime = time.time() - start_time

    outmessage = []
    outmessage.append("----------------------------------------------------------------------------------------")
    outmessage.append("ALLELIC ASSIGNMENT SUMMARY for \n%s"%(output_name_base))
    outmessage.append("----------------------------------------------------------------------------------------")
    outmessage.append("%d reads in Allele 1"%(a1_only))
    outmessage.append("%d reads in Allele 2"%(a2_only))
    outmessage.append("%d reads in Allele 1 with high probability"%(a1_count))
    outmessage.append("%d reads in Allele 2 with high probability"%(a2_count))
    outmessage.append("%d reads were not determined"%(nondetermined_count))
    outmessage.append("----- %s seconds -----" %(logtime))
    outmessage.append("%d BAD read names: "%(len(bad_reads)) + " , ".join(sorted(list(bad_reads))))
    outmessage.append("----------------------------------------------------------------------------------------")
    print("\n".join(outmessage))

    with open (logfile, "w") as f:
        f.write("\n".join(outmessage))


if __name__ == "__main__":
    #cProfile.run("main()")
    main()

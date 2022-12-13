#########################################################
#                                                       #
#               ASE Replicates Project                  #
#  https://github.com/gimelbrantlab/ASEReadCounter_star #
#                                                       #
#   Authors: Mendelevich Asia, Svetlana Vinogradova     #
#   Idea is taken from GATK pipeline (Stephane Castel)	#
#                                                       #
#########################################################
##
## NOTE: 
## 1. SAM-files should be already sorted by read names (samtools sort -n)
## 2. Read names must be non-emptyi
## 3. AS field must be present in SAM files (to filter by alignment quality) - typical for STAR
##
## USAGE: python3 alleleseparation.no_vA.py --samA1 [path to alignment on A1] --samA2 [path to alignment on A2] --odir [output dir] --obase [base name for output] --allele1 [name of the first allele] --allele2 [same for allele 2] --paired [0|1]
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

# special samtools ordering:
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

AS_prefix = "AS:i:"
def get_score(r_data, paired):
    if paired: astr = r_data[0]
    else: astr = r_data
    start_pos = astr.find(AS_prefix)+len(AS_prefix)
    fin_pos = astr.find('\t', start_pos)
    return int(astr[start_pos:fin_pos])

def get_name(line):
    return line[:line.find('\t')]

def get_name_and_score(line):
    name_pos = line.find("\t")
    score_start_pos = line.find("AS:i:", name_pos)+5
    score_fin_pos = line.find("\t", score_start_pos)
    return line[:name_pos], int(line[score_start_pos:score_fin_pos])

def asc_order(name1, name2):
    return (uporotiycmp(name1, name2) <= 0)

def main():
    # Parse arguments:

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

    # Open output_sam; Get header:

    header = subprocess.check_output("samtools view -SH "+ args.samA1, shell=True, universal_newlines=True)
    rg_a2  = subprocess.check_output("samtools view -SH "+ args.samA2 + '| grep "^@RG"', shell=True, universal_newlines=True)
    
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

    a1_only , a1_count , a2_only , a2_count , nondetermined_count = 0 , 0 , 0 , 0 , 0 
    bad_reads = set()

    source_A1 = open(args.samA1, 'r')
    source_A2 = open(args.samA2, 'r')

    # Skip header in each file:
    
    A1_skip = int(subprocess.check_output("samtools view -SH "+args.samA1+" | wc -l", shell=True, universal_newlines=True).strip())
    A2_skip = int(subprocess.check_output("samtools view -SH "+args.samA2+" | wc -l", shell=True, universal_newlines=True).strip())

    for i in range(m_skip):
        source_A1.readline()
    for i in range(p_skip):
        source_A2.readline()
    
    # Reads processing:

    def output_read(fhandler, a_read):
        if paired:
            fhandler.write(a_read[0])
            fhandler.write(a_read[1])
        else:
            fhandler.write(a_read)

    def blocks_generator(fhandler):
        beg_line = fhandler.readline()
        output = [beg_line]
        name, score = get_name_and_score(beg_line)
        our_prefix = name+'\t'
        for line in fhandler:
            if line.startswith(our_prefix):
                output.append(line)
            else:
                yield output, name, score
                beg_line = line
                output = [beg_line]
                name, score = get_name_and_score(beg_line)
                our_prefix = name+'\t'
        yield output, name, score

    def correct_blocks_generator(fhandler):
        if paired:
            for block, name, score in blocks_generator(fhandler):
                if (len(block) == 2): yield block, name, score
                else: bad_reads.add(name)
        else:
            for block, name, score in blocks_generator(fhandler):
                if (len(block) == 1): yield block[0], name, score
                else: bad_reads.add(name)
    
    # Create generator objects:
    
    A2gen = correct_blocks_generator(source_A2)
    A1gen = correct_blocks_generator(source_A1)
  
    # Separate till some EOF: 
    try:
        a1_read, a1_read_name, a1_score = A2gen.__next__()
        a2_read, a2_read_name, a2_score = A2gen.__next__()
        while 1:
            if a1_read_name == a2_read_name:
                if a1_score > a2_score:
                    output_read(out_stream_1a, a1_read)
                    a1_count += 1
                elif a1_score < a2_score:
                    output_read(out_stream_2a, a2_read)
                    a2_count += 1
                else:
                    nondetermined_count += 1
                    output_read(out_stream_Na, a1_read)
                a1_read, a1_read_name, a1_score = A1gen.__next__()
                a2_read, a2_read_name, a2_score = A2gen.__next__()
            elif asc_order(a2_read_name, a1_read_name):
                a2_only += 1
                output_read(out_stream_2a, a2_read)
                a2_read, a2_read_name, a2_score = A2gen.__next__()
            else:
                a1_only += 1
                output_read(out_stream_1a, a1_read)
                a1_read, a1_read_name, a1_score = A1gen.__next__()
    except StopIteration:
        pass

    # Write the remain part of reads:
    for a1_read, a1_read_name , a1_score in A1gen:
        a1_only += 1
        output_read(a1_read)
    for a2_read, a2_read_name , a2_score in A2gen:
        a2_only += 1
        output_read(a2_read)
            
    
    # Closing:
    out_stream_Na.close(); out_stream_1a.close(); out_stream_2a.close()
    source_A1.close(); source_A2.close()

    logfile = os.path.join(args.odir, ".".join([output_name_base,"log_allelic_assignment", "txt"]))
    logtime = time.time() - start_time

    # Writing the log:

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

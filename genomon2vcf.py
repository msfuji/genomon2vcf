#!/usr/bin/env python
"""
Convert SNV/INDEL calls by the Genomon2 pipeline into the VCF format.
The input is a single-sample variant call in mutation/ directory.
Merged calls in post_analysis/ shold not be used.
Usage: genomon2vcf.py genomon_mutation.result.filt.txt outfile [sample id]
"""
__author__ = "Masashi Fujita <m-fujita@riken.jp>"
__version__ = "0.0.1"
__date__    = "Sep. 1, 2018"

import sys
import datetime
import pysam

########################################################################

vcf_header="""\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
##fileDate=%s
##reference=%s
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s
""" 
vcf_header_date=datetime.datetime.now().strftime("%Y%m%d")

def die(msg):
    sys.stderr.write(msg+"\n")
    sys.exit(1)

genome_file="/home/w3varann/.genomon_local/genomon_pipeline-2.5.2/database/GRCh37/GRCh37.fa"

########################################################################

if len(sys.argv)==4:
    sample=sys.argv[-1]
    sys.argv.pop(-1)
else:
    sample="Tumor"

if len(sys.argv)!=3:
    die("Usage: genomon2vcf.py genomon_mutation.result.filt.txt outfile [sample id]")
infile=sys.argv[1]
outfile=sys.argv[2]

genome=pysam.FastaFile(genome_file)

with open(outfile, "w") as f_out:
    with open(infile) as f_in:
        # skip header
        for i in range(3):
            line=f_in.next()
            if line[0] != "#":
                die("ERROR: invalid file format\n")
        line=f_in.next()
        if line.split()[0] != "Chr":
            die("ERROR: invalid file format\n")

        f_out.write(vcf_header % (vcf_header_date, genome_file, sample))
        for line in f_in:
            a=line.strip("\n").split("\t")
            CHROM=a[0]
            POS=a[1]
            ID="."
            REF=a[3]
            ALT=a[4]
            QUAL="."
            FILTER="PASS"
            INFO="SOMATIC"
            FORMAT="GT"
            GT="0/1"

            ########## Insertion ##########
            if REF=="-":
                p=int(POS)
                base=genome.fetch(CHROM, p-1, p)
                REF=base
                ALT=base+ALT
            ########## Deletion ##########
            elif ALT=="-":
                p=int(POS) - 1
                POS=str(p)
                base=genome.fetch(CHROM, p-1, p)
                REF=base+REF
                ALT=base

            newline="\t".join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, GT])
            f_out.write(newline+"\n")
        

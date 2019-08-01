#!/usr/bin/env python
"""
Convert SNV/INDEL calls by the Genomon2 pipeline into the VCF format.
The input is a single-sample variant call in mutation/ directory.
Merged calls in post_analysis/ should not be used.
Usage: python genomon2vcf.py genomon_mutation.result.filt.txt
"""
__author__ = "Masashi Fujita <mssfjt@gmail.com>"
__version__ = "0.0.1"
__date__ = "July 31, 2019"

import sys
import argparse
import pysam
import pandas as pd


def correct_position(row, genome):
    #
    # Insertion
    #
    if row.Ref == '-':
        vcf_pos = row.Start
        base = genome.fetch(row.Chr, vcf_pos - 1, vcf_pos)
        vcf_ref = base
        vcf_alt = base + row.Alt
    #
    # Deletion
    #
    elif row.Alt == '-':
        vcf_pos = row.Start - 1
        base = genome.fetch(row.Chr, vcf_pos - 1, vcf_pos)
        vcf_ref = base + row.Ref
        vcf_alt = base
    #
    # SNV
    #
    else:
        vcf_pos = row.Start
        vcf_ref = row.Ref
        vcf_alt = row.Alt

    return pd.Series({'POS': vcf_pos, 'REF': vcf_ref, 'ALT': vcf_alt})


def build_info(row):
    return "VAF={:.3f}".format(row.misRate_tumor)


#
# parse args
#
parser = argparse.ArgumentParser(description='Convert Genomon2 SNV/INDEL calls to VCF.')
parser.add_argument('infile', metavar='genomon_file', help='Genomon2 SNV/INDEL file for a single sample')
parser.add_argument('--out', '-o', metavar='VCF', help='Output VCF file [default: stdout]')
parser.add_argument('--ref', '-r', metavar='FAS', help='FASTA file of reference human genome')
parser.add_argument('--sample', '-s', metavar='ID', help='Sample ID')
args = parser.parse_args()

#
# prepare VCF header
#
vcf_header="""\
##fileformat=VCFv4.1
##FORMAT=<ID=VAF,Number=1,Type=Fload,Description="Variant allele frequency">
"""
if args.ref is not None:
    vcf_header += '##reference="' + args.ref + '"\n'
vcf_header += """#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""
if args.sample is not None:
    vcf_header += '\t' + args.sample
vcf_header += '\n'

#
# open input
#
df = pd.read_csv(args.infile, sep='\t', skiprows=3)



#
# Correct INDEL positions if ref genome is given
#
if args.ref is not None:
    genome = pysam.FastaFile(args.ref)
    corrected_pos = df.apply(lambda row: correct_position(row, genome), axis=1)
    df = pd.concat([df, corrected_pos], axis=1)
else:
    df['POS'] = df.Start
    df['REF'] = df.Ref
    df['ALT'] = df.Alt

#
# format into VCF
#
df['ID'] = '.'
df['QUAL'] = '.'
df['FILTER'] = '.'
df['INFO'] = df.apply(build_info, axis=1)
df = df[['Chr', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]

#
# save
#
if args.out is not None:
    f_out = open(args.out, 'w')
else:
    f_out = sys.stdout
f_out.write(vcf_header)  # write VCF header
df.to_csv(f_out, sep='\t', header=False, index=False)
f_out.close()
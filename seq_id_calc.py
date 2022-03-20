'''
04/11/2019
Python script to calculate the sequence identity score for CIGAR and FATCIGAR strings for each read in a BAM file. 
'''

import sys
import argparse
import pysam
import decimal
try:
    from itertools import groupby, imap
except ImportError:
    imap=map
import re


def cigar_seqid(cigarstring):
    """
    Given the CIGAR string, it returns the sequence identity score based on the number of matched bases (includes mismatched bases) over the total number of bases.
    :param str cigarstring: string to calculate the identity scores from
    :rtype: decimal
    """
    #substitutes out any deletions and hard-clips and keeps just the digits from the CIGAR string
    cig_split = map(int, ["".join(x) for is_number, x in groupby(re.sub(r'\d+[D|H]', '', cigarstring), key=str.isdigit) if is_number is True])
    #totals up all the values
    total = sum(cig_split)
    #extracts just the number of matches from the CIGAR string
    match_split = map(int, ["".join(x) for is_number, x in groupby(re.sub(r'\d+[D|H|I|S]', '', cigarstring), key=str.isdigit) if is_number is True])
    #calculates the total number of matches
    match = sum(match_split)
    #returns the sequence identity score rounded up to 3 d.p
    return (decimal.Decimal(match)/decimal.Decimal(total)).quantize(decimal.Decimal("1.000"))


def fatcigar_seqid(fatcigarstring):
    """
    Given the FATCIGAR string, it returns the sequence identity score based on the number of truly matched bases over the total number of bases.
    :param str fatcigarstring: string to calculate the identity scores from
    :rtype: decimal
    """
    #substitutes out any deletions and hard-clips and keeps just the digits from the FATCIGAR string
    cig_split = map(int, ["".join(x) for is_number, x in groupby(re.sub(r'\d+[D|H]', '', fatcigarstring), key=str.isdigit) if is_number is True])
    #totals up all the values
    total = sum(cig_split)
    #extracts just the number of matches from the FATCIGAR string
    match_split = map(int, ["".join(x) for is_number, x in groupby(re.sub(r'\d+[D|H|I|S|X]', '', fatcigarstring), key=str.isdigit) if is_number is True])
    #calculates the total number of matches
    match = sum(match_split)
    #returns the sequence identity score rounded up to 3 d.p
    return (decimal.Decimal(match)/decimal.Decimal(total)).quantize(decimal.Decimal("1.000"))


def main():
    #Instantiate the parser
    parser = argparse.ArgumentParser(description='Calculates sequence identity scores from the CIGAR and FATCIGAR string for sequence reads in a BAM file.')

    #Required positional argument
    parser.add_argument('bam_file', type=str,
                    help='A required string argument stating the name of the input BAM file')
    parser.add_argument('outbam_file', type=str,
                    help='Output BAM file to write the reads with the sequence identity scores. The CIGAR string scores will be written to the ZA tag and FATCIGAR string to the ZB tag')

    #Parse arguments
    args = parser.parse_args()

    #Read in BAM files and output file
    bam = pysam.AlignmentFile(args.bam_file, "rb")
    outbam = pysam.AlignmentFile(args.outbam_file, "wb", template=bam)

    #reads in the lines from each BAM file
    for read in bam.fetch():
        if not read.is_unmapped:
            #if read is mapped, calculates the sequence identity score for both the CIGAR and FATCIGAR string and writes the scores to the output file in tab separated format
            read.set_tag("ZA", int(cigar_seqid(read.cigarstring)))
            read.set_tag("ZB", int(fatcigar_seqid(read.get_tag("XG"))))
        outbam.write(read)


    #close files
    bam.close()
    outbam.close()

if __name__== "__main__":
  main()

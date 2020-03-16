#!/usr/bin/env python3
import sys
import pysam

def main(argv):

    # parameters
    filename = argv[1]
    positions = [int(item) for item in argv[2:]]

    # Open valid input files 
    if '.sam' in filename.lower():
        samfile = pysam.AlignmentFile(filename, "r")
    elif '.bam' in filename.lower():
        pysam.index(filename)
        samfile = pysam.AlignmentFile(filename, "rb")
    else:
        print('Invalid input file.')
        return

    # Read data
    lines = samfile.fetch()

    align = [line for line in lines]
    align = [str(line).strip().split('\t') for line in align]
    align = [[int(line[3]), line[5], line[9]] for line in align] # pos, cigar, query

    samfile.close()  

if __name__ == '__main__':
    main(sys.argv)
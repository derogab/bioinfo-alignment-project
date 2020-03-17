#!/usr/bin/env python3
import sys
import pysam
import argparse

def output(data):
    print('OUTPUT DATA: \n')

    if len(data) == 0:
        print('No reads found.')
    else:
        for x in data:
            print(str(x[0]) + '\t' + x[2] + '\n') # pos, query

def main(argv):

    # Program info and argument parser
    parser = argparse.ArgumentParser(
        description='Print reads from a SAM/BAM file filtered for start positions.'
    )

    parser.add_argument(
        '-f', '--file',
        required=True,
        help='path to input file'
    )

    parser.add_argument(
        '-p', '--pos',
        nargs='+',
        type=int,
        help='list of positions'
    )

    args = parser.parse_args()

    # parameters
    filename = args.file
    positions = args.pos

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

    align = [str(line).strip().split('\t') for line in lines]
    align = [[int(line[3]), line[5], line[9]] for line in align] # pos, cigar, query

    # Filter
    if args.pos is None:
        align = [line for line in align if line[0] in positions]
    
    # Output requested data 
    output(align)

    # Close file 
    samfile.close()  

if __name__ == '__main__':
    main(sys.argv)
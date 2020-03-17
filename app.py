#!/usr/bin/env python3
import sys
import pysam
import argparse

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

    # Open valid input files 
    if '.sam' in args.file.lower():
        samfile = pysam.AlignmentFile(args.file, "r")
    elif '.bam' in args.file.lower():
        pysam.index(args.file)
        samfile = pysam.AlignmentFile(args.file, "rb")
    else:
        print('Invalid input file.')
        return

    # Read data
    lines = samfile.fetch()

    align = [str(line).strip().split('\t') for line in lines]

    # Filter
    if args.pos is not None:
        align = [line for line in align if int(line[3]) in args.pos]
    
    # Output requested data 
    print('OUTPUT DATA: \n')

    if len(align) == 0:
        print('No reads found.')
    else:
        for x in align:
            print(x[9] + ' (' + x[3] + ') \n') # query (pos)

    # Close file 
    samfile.close()  

if __name__ == '__main__':
    main(sys.argv)
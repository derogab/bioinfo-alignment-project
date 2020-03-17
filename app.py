#!/usr/bin/env python3
import sys
import pysam

def output(data):
    print('OUTPUT DATA: \n')

    if len(data) == 0:
        print('No reads found.')
    else:
        for x in data:
            print(str(x[0]) + '\t' + x[2] + '\n') # pos, query

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

    align = [str(line).strip().split('\t') for line in lines]
    align = [[int(line[3]), line[5], line[9]] for line in align] # pos, cigar, query

    # Filter 
    align = [line for line in align if line[0] in positions]
    
    # Output requested data 
    output(align)

    samfile.close()  

if __name__ == '__main__':
    main(sys.argv)
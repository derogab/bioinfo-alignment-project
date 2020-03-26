#!/usr/bin/env python3
import sys
import pysam
import argparse
import numpy as np

# Global variables
result_counter = 0

# Get positions range length
def get_pos_range_length(line):
    seq_length = 0
    for op in line.cigar:
        seq_length += op[1]
    return seq_length


# Output template
def print_read(line):
    global result_counter

    # print a line / template output 
    print(line.query_sequence + ' (' + str(line.reference_start) + ' - ' + str(line.cigarstring) + ' - ' + str(get_pos_range_length(line)) + ')')

    # increment counter value
    result_counter += 1


# Main function w/ arguments 
def main(argv): 
    global result_counter

    # Program info
    parser = argparse.ArgumentParser(
        description='Print reads from a SAM/BAM file filtered for start positions.'
    )
    
    # Main features arguments 
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

    # Parse arguments
    args = parser.parse_args()

    # Open valid input files 
    if '.sam' in args.file.lower():
        open_mode = 'r'
    elif '.bam' in args.file.lower():
        open_mode = 'rb'
        # create index file for binary file
        pysam.index(args.file)
    else:
        print('Invalid input file.')
        return

    try:
        samfile = pysam.AlignmentFile(args.file, open_mode)
    except Exception as ex:
        print('An exception of type ' + type(ex).__name__ + ' occurred: ' + ex.args[1])
        return

    # Core
    if args.pos is None:
        
        # get all queries 
        for line in samfile.fetch():
            print_read(line)

    else: 

        # get only queries of selected positions 
        for line in samfile.fetch():
            if(np.intersect1d(list(range(line.pos, line.pos + get_pos_range_length(line) - 1)), args.pos)):
                print_read(line)

    # Counter results
    if result_counter == 0:
        print('No results found.')
    else:
        print(str(result_counter) + ' results printed.')

    # Close file 
    samfile.close()


if __name__ == '__main__':
    main(sys.argv)
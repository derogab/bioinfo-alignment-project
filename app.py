#!/usr/bin/env python3
import sys
import pysam
import argparse

# Global variables
result_counter = 0

# Functions 
def range_control(items, limits):
    # Check if almost one item is between limits
    for item in items:
        if item in limits:
            return True
    return False

def print_read(line):
    global result_counter

    # Output template 
    print('len(seq):\t', line.query_length)
    print('pos:\t\t', line.reference_start)
    print('cigar:\t\t', line.cigarstring)
    print('subreference:\t', line.query_sequence)
    print('query:\t\t', line.query_alignment_sequence)
    print()

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

    # Core feature 
    if args.pos is None:
        
        # get all queries 
        for line in samfile.fetch():
            print_read(line)

    else: 

        # get only queries of selected positions 
        for line in samfile.fetch():
            if(range_control(args.pos, range(line.pos, line.pos + line.query_length))):
                print_read(line)

    # Counter results
    if result_counter == 0:
        print('No results found.')
    else:
        print(str(result_counter) + ' results printed.')

    # Close file 
    samfile.close()

# Start main  w/ arguments 
if __name__ == '__main__':
    main(sys.argv)
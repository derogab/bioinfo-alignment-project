#!/usr/bin/env python3
import sys
import pysam
import argparse

# Global variables
result_counter = 0

# Output template
def print_read(line):
    global result_counter
    
    # print a line / template output 
    print(line.query_alignment_sequence + ' (' + str(line.reference_start) + ')')

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

    # Extra features arguments
    parser.add_argument(
        '--only-pos',
        action='store_true',
        help='get list of all align start positions'
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


    # Main features
    if args.pos is None and not args.only_pos:
        
        # get all queries 
        for line in samfile.fetch():
            print_read(line)

    elif args.pos is not None and not args.only_pos: 

        # get only queries of selected positions 
        for line in samfile.fetch():
            if line.reference_start in args.pos:
                print_read(line)

    # Extra features
    elif args.pos is None and args.only_pos:

        # get all positions
        for line in samfile.fetch():
            print(line.reference_start)
            result_counter += 1

    elif args.pos is not None and args.only_pos: 

        # get only selected positions
        for line in samfile.fetch():
            if line.reference_start in args.pos:
                print(line.reference_start)
                result_counter += 1

    # Counter results
    if result_counter == 0:
        print('No results found.')
    else:
        print(str(result_counter) + ' results printed.')

    # Close file 
    samfile.close()


if __name__ == '__main__':
    main(sys.argv)
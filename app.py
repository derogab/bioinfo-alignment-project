#!/usr/bin/env python3
import sys
import pysam
import argparse

def main(argv):

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
        output = [line.query_alignment_sequence for line in samfile.fetch()]
    
    elif args.pos is not None and not args.only_pos: 

        # get only queries of selected positions 
        output = [line.query_alignment_sequence for line in samfile.fetch() if line.reference_start in args.pos]

    # Extra features
    elif args.pos is None and args.only_pos:

        # get all positions 
        output = [line.reference_start for line in samfile.fetch()]

    elif args.pos is not None and args.only_pos: 

        # get only selected positions
        output = [line.reference_start for line in samfile.fetch() if line.reference_start in args.pos]


    # Output requested data 
    print('OUTPUT DATA: \n')

    if len(output) == 0:
        print('No reads found.')
    else:
        print(output)

    # Close file 
    samfile.close()


if __name__ == '__main__':
    main(sys.argv)
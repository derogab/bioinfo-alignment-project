#!/usr/bin/env python3
import sys
import pysam
import argparse

def main(argv):

    # Program info
    parser = argparse.ArgumentParser(
        description='Print reads from a SAM/BAM file filtered for start positions.'
    )
    
    # Main features
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

    # Extra features
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

    # Read data
    lines = samfile.fetch()

    align = [str(line).strip().split('\t') for line in lines]

    # Get only positions available
    if(args.only_pos):
        positions_list = [int(line[3]) for line in align]
        print(positions_list)
        return

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
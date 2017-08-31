#! /usr/bin/env python 
"""
Read in a file with a record of the chemical abundance
for each spectral line, and comment out the
rows with non-numerical data.
"""

import sys
import argparse

def comment_the_file(input_filepath, output_path):
    """Comments rows in the file with non-numerical data."""

    f_in = open(input_filepath,'r')
    data = f_in.read()
    lines = data.split('\n')
    f_in.close()

    out_lines = []
    for line in lines:
	columns = line.split()
        try:
           col1 = float(columns[0])
	   col2 = float(columns[1])
	   col3 = float(columns[2])
	   col4 = float(columns[3])
	   col5 = float(columns[4])
	   col6 = float(columns[5])
	   col7 = float(columns[6])
           line_format = '%10.2f %7.2f %10.3f %7.1f %8.2f %7.2f %8.2f'
           out_lines.append(line_format % (col1,
                                           col2,
                                           col3,
                                           col4,
                                           col5,
                                           col6,
                                           col7))
	except IndexError:
	   continue 
	except ValueError:
	   out_lines.append('# '+line)

    outtext = '\n'.join(out_lines)
    if output_path == sys.stdout:
       # Add a newline character to the end of the outtext
       # string to make the terminal output cleaner.
       output_path.write(outtext+'\n')
    else:
       with open(output_path, 'w') as f_out:
            f_out.write(outtext)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_filepath',
                         type=str,
                         help='Path of the input file.')
    parser.add_argument('-o', '--output_path',
                         type=str,
                         default=sys.stdout,
                         help='Path for the destination file. If no path '
                              'is specified, output is written to stdout. '
                              'Note that if a file exists at this path, '
                              'it will be overwritten.')
    args = parser.parse_args()
    comment_the_file(args.input_filepath, args.output_path)

#!/usr/bin/env python
"""
Convert two data files to a format that can be inserted into a 
LaTeX table. 
"""

import numpy as np
import sys
import argparse


def get_elem_name(elem_number):
    """Converts a floating point number representing an element to a string name.
      
    The integer component of the number should be the element's atomic number,
    and the decimal component should be its ionized state (i.e, 1 = singly ionized,
    2 = doubly ionized, etc.)
    """

    # Currently supports only neutral and singly ionized iron
    elem_name='Fe'
    if (elem_number != 26.1):
        elem_ion='I'
    else:
        elem_ion='II'
    return elem_name, elem_ion


def reformat_files(datafile1_path, datafile2_path, output_path):
    """Converts data from two files into a LaTeX table format.

    The data files should be modified from the raw output, i.e.,
    all rows with non-numerical data should be commented with a '#'.
    For the variable names:
         wave = wavelength
         elem = floating point numerical representation of the element
         ep = excitation potential for the transition
         loggf = log of the oscillator strength for the transition
         ew = equivalent width of the spectral line
         logN = log of the abundance of the element
    """

    data1_wave, data1_elem, data1_ep, data1_loggf, data1_ew, data1_logN = \
               np.loadtxt(datafile1_path, 
                          usecols=(0,1,2,3,4,6), 
                          unpack=True)

    # For datafile2, the values for wave, elem, ep, and loggf should
    # be identical.
    data2_ew, data2_logN = \
               np.loadtxt(datafile2_path,
                          usecols=(4,6),
                          unpack=True)

    latex_table_format = '\\ion{%2s}{%1s} & %8.2f & %5.2f & %7.3f & %5.1f & %5.2f & %5.1f & %5.2f \\\\'

    out_lines = []
    for ix, wave_item in enumerate(data1_wave):
        elem_name, elem_ion = get_elem_name(data1_elem[ix])
        out_lines.append(latex_table_format % (elem_name,
                                               elem_ion,
                                               wave_item,
                                               data1_ep[ix],
                                               data1_loggf[ix],
                                               data1_ew[ix],
                                               data1_logN[ix],
                                               data2_ew[ix],
                                               data2_logN[ix]))
    outtext = '\n'.join(out_lines)
    if output_path == sys.stdout:
       # Add a newline character to the end of the outtext
       # string to make the terminal output cleaner.
       output_path.write(outtext+'\n')
    else:
       with open(output_path, 'w') as f:
	    f.write(outtext)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('data1_filepath',
                         type=str,
                         help='Path of the input file for the 1st star.')
    parser.add_argument('data2_filepath',
                         type=str,
                         help='Path of the input file for the 2nd star '
                              '(usu. a reference star like the Sun).')
    parser.add_argument('-o', '--output_path',
                         type=str,
                         default=sys.stdout,
                         help='Path for the destination file. If no path '
                              'is specified, output is written to stdout. '
                              'Note that if a file exists at this path, '
                              'it will be overwritten.')
    args = parser.parse_args()
    reformat_files(args.data1_filepath, args.data2_filepath, args.output_path)

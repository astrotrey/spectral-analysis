#! /usr/bin/env python 
#
# NAME:
#   modify_abd.py
#
# PURPOSE:
#   Read in a file with a record of the chemical abundance
#   for each spectral line, and comment out the
#   rows with non-numerical data.
#
# CALLING SEQUENCE:
#
# ./modify_abd.py input_filepath output_filepath
#
#
import math
import sys

input_filepath = sys.argv[1]
output_filepath = sys.argv[2]


fid = open(input_filepath,'r')
data = fid.read()
lines = data.split('\n')
fid.close()


#-Add a '#' to the start of each
#-line of text in the input file
#-
for line in lines:

	columns = line.split()

        try:
           col1 = float(columns[1])
	   col2 = float(columns[2])
	   col3 = float(columns[3])
	   col4 = float(columns[4])
	   col5 = float(columns[5])
	   col6 = float(columns[6])
	   col7 = float(columns[7])

           fid = open(output_filepath,'a')
           form='%10.2f %7.2f %10.3f %7.1f %8.2f %7.2f %8.2f \n'
	   fid.write(form % (col1,col2,col3,col4,col5,col6,col7))
           fid.close()
	   
	except IndexError:
	   continue 
	except ValueError:
           fid = open(output_filepath,'a')
	   fid.write('# '+line+'\n')
           fid.close()

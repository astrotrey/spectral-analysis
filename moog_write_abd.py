#! /usr/bin/env python 
#
# NAME:
#   moog_write_abd.py
#
# PURPOSE:
#   Read in a moog output file, and comment out the
#   lines with non-numerical data.
#
# CALLING SEQUENCE:
#
# ./write_moog_abd.py moog_abd_file outfile
#
#
import math
import sys

file = sys.argv[1]
outfile = sys.argv[2]


#-Add a '#' to the start of each
#-line of text in the input file
#-
fid = open(file,'r')
data = fid.read()
lines = data.split('\n')

fid = open(outfile,'w')
form = '%1s\n'
fid.write(form % ('#'))

for line in lines:

        fid = open(outfile,'a')
        fid.write('# '+line+'\n')



#-Strip off the '#' for 
#-the lines with numerical data
#-
fid = open(outfile,'r')
data = fid.read()
lines = data.split('\n')

fid = open(outfile,'w')
form = '%1s\n'
fid.write(form % ('#'))

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

           fid = open(outfile,'a')
           form='%10.2f %7.2f %10.3f %7.1f %8.2f %7.2f %8.2f \n'
	   fid.write(form % (col1,col2,col3,col4,col5,col6,col7))
	   
	except IndexError:
	   continue 
	except ValueError:
           fid = open(outfile,'a')
	   fid.write(line+'\n')


#END OF PROGRAM

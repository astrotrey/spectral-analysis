#!/usr/bin/env python
#
#
#
#
"""
Convert moog abd output files to latex-formatted
table rows that can be copied and pasted into latex
table environment
"""


#--###################################################
#-- Import modules and parse the input arguments
#--
#--###################################################
import numpy as np
import sys


#-- Print the command and arguments to the screen
#--
print(sys.argv)

#-- Assign the input arguments to variables
#--
ifile01 = sys.argv[1]
ifile02 = sys.argv[2]
ofile01 = sys.argv[3]


#--###################################################
#-- Read in the data from the (modified) moog abd files
#--
#--###################################################


#-- Read in the data from first input file
#--
wave01, elem01, ep01, loggf01, ew01, logN01 = np.loadtxt(ifile01, usecols=(0,1,2,3,4,6), unpack=True)
#print wave01, elem01

#-- Read in the wavelength from second input file as way to check that the
#-- the data arrays in the two files have the same dimension
#--
wave02, ew02, logN02 = np.loadtxt(ifile02, usecols=(0,4,6), unpack=True)
print wave02, ew02, logN02



#--###################################################
#-- Write the data to ofile01 in latex-table format
#--
#--###################################################

#-- Initialize variables for the loop,
#-- and open the output file for writing,
#-- so that any pre-existing version of the
#-- file is deleted
#--
i=0
fid=open(ofile01,'w')

while (i <= len(wave01)-1):

    #-- Determine the text string for
    #-- the variable elem_name given the
    #-- numerical value of elem01[i]
    #-- (Next: try to use a dictionary,
    #--        case statement, or sequence
    #--        of elif statements to accomplish
    #--        this in a more general way for all
    #--        chemical species. Elif or case statement
    #--        implementation would probably be better
    #--        handled by a pre-defined function.)
    #--
    if (elem01[i] !=  26.1):
        elem_name='Fe'
        elem_ion='I'
    else:
        elem_name='Fe'
        elem_ion='II'


    #-- Define the format and print to the
    #-- output text file row-by-row
    #-
    fid=open(ofile01,'a')
    form='\\ion{%2s}{%1s} & %8.2f & %5.2f & %7.3f & %5.1f & %5.2f & %5.1f & %5.2f \\\\ \n'
    fid.write(form % (elem_name,elem_ion,wave01[i],ep01[i],loggf01[i],ew01[i],logN01[i],ew02[i],logN02[i]))

    i=i+1





### END OF PROGRAM

#! /usr/bin/env python 
"""
Evaluate the chemical model for a
given star to see how well it fits
the data.
"""

import math
import cmath
import numpy as np
import sys
import argparse
import subprocess
from abd2latex import get_elem_name

def get_params(param_filepath):
    """Extract the model's parameters from the parameter file."""
    #return (Teff, logg, feh, vmicro, atm_count, logfile, logfile2, logfile3,
    #        datafile(?), model_filepath)

def read_datafile(data_filepath):
    """Extract the relevant data from a given star file.

       Note on variable names:
            wave = wavelength
            ep = excitation potential
            ew = equivalent width
            rew = reduced equivalent width
            logN = the absolute chemical abundance
    """

    f_in = open(data_filepath, 'r')
    data = f_in.read()
    lines = data.split('\n')
    f_in.close()

    wave = []
    elem_num = []
    ep = []
    ew = []
    rew = []
    logN = []
    for line in lines:
        columns = line.split()
        try:
           wave.append(float(columns[0]))
           elem_num.append(float(columns[1]))
           ep.append(float(columns[2]))
           ew.append(float(columns[4]))
           rew.append(float(columns[5]))
           logN.append(float(columns[6]))
        except IndexError:
	   continue
	except ValueError:
           continue
  
    return wave, elem_num, ep, ew, rew, logN

def check_model(data1_filepath, data2_filepath, param_filepath):
    """Determine if the model is a good fit to the data."""
     
    star1_wave, star1_elem_num, star1_ep, star1_ew, star1_rew, star1_logN =\
                read_datafile(data1_filepath)
    star2_wave, star2_elem_num, star2_ep, star2_ew, star2_rew, star2_logN =\
                read_datafile(data2_filepath)
    
    if star1_wave != star2_wave:
       print '\nERROR! The spectral line data from the two input files are not identical.\n'
       print 'Check '+data1_filepath+' to '+data2_filepath+' make sure they are the correct files.\n'
       exit()

    star1_fe1_ep = []   # ep, ew, and rew needed for statistics that
                        # determine how well the model fits the data.
                        # But only needed for star1 and Fe-I.
    star1_fe1_ew = []
    star1_fe1_rew = []
    star1_fe1_logN = []
    star1_fe2_logN = [] # Only need the absolute abundances of Fe-II
                        # from the 1st star.
    star2_fe1_logN = [] # 2nd star is reference; only need Fe-I and Fe-II logNs
    star2_fe2_logN = [] # So Yeah: why read it in if you don't need it? Edit read_datafile()??
    for ix, wave_item in enumerate(star1_wave):
        #elem_name, elem_ion = get_elem_name(star1_elem_num[ix]) #might use this in a later version
        if (star1_elem_num[ix] != 26.1):
            star1_fe1_ep.append(star1_ep[ix])
            star1_fe1_ew.append(star1_ew[ix])
            star1_fe1_rew.append(star1_rew[ix])
            star1_fe1_logN.append(star1_logN[ix])
            star2_fe1_logN.append(star2_logN[ix])
        else:
            star1_fe2_logN.append(star1_logN[ix])
            star2_fe2_logN.append(star2_logN[ix])

    # Evaluate how well the model fits the data.
    # The model is considered a solution if it fits the data well enough,
    # such that:
    #      1) star1_fe1_avg - star1_fe2_avg = 0.00
    #             The Fe-I and Fe-II abundances should agree when rounded
    #             to the nearest hundredth.
    #      2) fe1_vs_ep < 0.01 
    #             The correlation coefficient for star1_fe1 vs star1_fe1_ep
    #             must be less than 0.01 (no rounding!)
    #      3) fe1_vs_rew < 0.01
    #             The correlation coefficient for star1_fe1 vs star1_fe1_rew
    #             must be less than 0.01 (no rounding!)
    #
    star1_fe1 = np.subtract(star1_fe1_logN, star2_fe1_logN)  # Fe-I relative abundances of star1
                                                             # relative to star2
    star1_fe1_avg = np.mean(star1_fe1)                       # Fe-I mean relative abundance
    print star1_fe1
    print star1_fe1_avg

    star1_fe2 = np.subtract(star1_fe2_logN, star2_fe2_logN)  # Fe-II relative abundances of star1
                                                             # relative to star2
    star1_fe2_avg = np.mean(star1_fe2)                       # Fe-II mean relative abundance
    print star1_fe2
    print star1_fe2_avg
        
    fe1_vs_ep = np.corrcoef(star1_fe1, star1_fe1_ep)         #corr. coeff. for Fe-I vs ep
    fe1_vs_rew = np.corrcoef(star1_fe1, star1_fe1_rew)       #corr. coeff. for Fe-I vs rew
    print fe1_vs_ep 
    print fe1_vs_rew
    

    # Is subprocess really needed to delete the files when the model is not a good
    # fit to the data?


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('data1_filepath',
                         type=str,
                         help='Path of data file for the 1st star.')
    parser.add_argument('data2_filepath',
                         type=str,
                         help='Path of the data file for the 2nd star '
                              '(usu. a refence star like the Sun).')
    parser.add_argument('param_filepath',
                         type=str,
                         help='Path of the param file for the model: '
                              'must contain Teff, logg, etc for the model.')
    args = parser.parse_args()
    check_model(args.data1_filepath, args.data2_filepath, args.param_filepath)

################################################################################################
################################################################################################
################################################################################################

"""
#- Read value for ew_corr option
#-
ew_corr = int(sys.argv[3])

#- Read in Teff, log_g, fe_h, and vmicro
#- from the model atmosphere
#-
n_Fe1 = float(sys.argv[11])
n_Fe2 = float(sys.argv[12])
star_file_mod = (sys.argv[13]) #starfile deleted if parameters are not a solution
star_log = (sys.argv[14]) #starlog deleted if parameters are not a solution
star_model = (sys.argv[15]) #starmodel deleted if parameters are not a solution
logfile3 = 'soluxns.log'
logfile4 = 'soluxns_almost.log'
fid1 = open(logfile2,'w')

#- Also count the number of lines
#- measured for each element
#-
elem_strs = [0,0]
elem_n_lines = [n_Fe1,n_Fe2]


              if (logfile2 != 'foo.log'):
	          form='%7.2f  %9.2f %7.2f %7.2f %7.2f\n' 
	          fid1.write(form % (lambda1,lambda2,logN1,logN2,star_Fe2))
      
if (logfile2 != 'foo.log'):
   fid1.write(' \n')
   fid1.write('[FeI/H] = %6.3f\n' % (star_Fe1_avg))
   fid1.write('[FeII/H] = %6.3f\n' % (star_Fe2_avg))
   form='|[FeI/H]-[FeII/H]| = %6.3f\n'
   fid1.write(form % (diff))

#-Check if abs([Fe1/H] - [Fe2/H]) < 0.005
#-when [Fe1/H] and [FeII/H] are rounded to
#-the nearest hundredth
#-
digit3_Fe1 = star_Fe1_avg*10**(3) % 10
digit3_Fe2 = star_Fe2_avg*10**(3) % 10
Fe1 = star_Fe1_avg
Fe2 = star_Fe2_avg

#-If the value of the third decimal place is greater
#-than or equal to five, then round up to the nearest hundredth
#-(e.g., 0.345 -> 0.35)
#-this way, (Fe1=0.345,Fe2=0.349) is a solution but
#-(Fe1=0.344,Fe2=0.349) is not
#-
if (digit3_Fe1 >= 5):
	Fe1 = star_Fe1_avg - digit3_Fe1*10**(-3) + 10**(-2)

if (digit3_Fe2 >= 5):
	Fe2 = star_Fe2_avg - digit3_Fe2*10**(-3) + 10**(-2)

diff_Fe1Fe2 = Fe1 - Fe2

#-Determine the correlations for
#-[FeI/H] vs EP and [FeI/H] vs REW
#-
#-Also do EW vs EP if ew_corr==1
#
if ew_corr == 1:
 
   fid1.write(' \n')
   form = 'Pearson corr. coeff. of EW vs EP for Fe1 = %.4f\n'
   fid1.write(form % (r_EWvsEP_Fe1))

outfile=logfile
fid2=open(outfile,'a')
form='%d %6.2f %6.2f %6.2f %6.3f %6.3f %9.4f %9.4f \n'
fid2.write(form % (Teff,logg,feh,vmicro,star_Fe1_avg,star_Fe2_avg,r_Fe1vsEP,r_Fe1vsREW))


#soluxn_count=0
if (abs(diff_Fe1Fe2) < 0.005) and  (abs(r_Fe1vsEP) < 0.01) and  (abs(r_Fe1vsREW) < 0.01):
        #soluxn_count=soluxn_count+1
	print "\n\n YOU'VE FOUND A SOLUTION!\n"
	fid2.write('The above set of stellar parameters is a solution!\n')
	#sys.exit(1)
        args = ['touch', 'soluxn_found.tmp']
        subprocess.Popen(args)

        ofile=logfile3
        fid3=open(ofile,'a')
        fid3.write(form % (Teff,logg,feh,vmicro,star_Fe1_avg,star_Fe2_avg,r_Fe1vsEP,r_Fe1vsREW))

        args = ['rm', '-f', star_file] #-deleting *Fe_mod.abd
        subprocess.Popen(args)
        #-NOTE: remember that star_file is *Fe_mod.abd and star_file_mod is just *Fe.abd
        #-      maybe this should be changed in the future.

elif (abs(diff_Fe1Fe2) < 0.02) and  (abs(r_Fe1vsEP) < 0.01) and  (abs(r_Fe1vsREW) < 0.01):
        #-write the stellar parameters that are almost a solution to logfile4
        ofile=logfile4
        fid4=open(ofile,'a')
        fid4.write(form % (Teff,logg,feh,vmicro,star_Fe1_avg,star_Fe2_avg,r_Fe1vsEP,r_Fe1vsREW))

	fid2.write('The above set of stellar parameters is almost a solution.\n')
        args = ['touch', 'soluxn_almost_found.tmp']
        subprocess.Popen(args)

        args = ['rm', '-f', star_file, star_file_mod, star_log, star_model]
        subprocess.Popen(args)
        #-deleting '*.abd,*.log,*.atm' for the given set of stellar parameters

else:
        args = ['rm', '-f', star_file, star_file_mod, star_log, star_model]
        subprocess.Popen(args)
        #-deleting '*.abd,*.log,*.atm' for the given set of stellar parameters
        
"""

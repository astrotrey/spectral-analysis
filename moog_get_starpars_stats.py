#! /usr/bin/env python 
#
# NAME:
#
#   moog_get_starpars_stats.py
#
#
# PURPOSE:
#
#   Determine [FeI/H] and [FeII/H]
#   and Pearson correlation for
#   [FeI/H] vs EP and [FeI/H] vs REW
#
#
# CALLING SEQUENCE:
#
#   ./moog_get_starpars.py star_file sun_file
#
#
# INPUTS:
#
#   star_file - moog.abd file for a given
#               model atm 
#
#   sun_file - moog.abd file for the sun
#              (Teff=5777,logg=4.44,[Fe/H]=0.00,microt=1.38)
#
#   ew_corr - determines if r_EWvsEP is printed to output   
#
#
# Last revised: 2015/09/25 by Trey Mack
#
#

import math
import cmath
import sys
import subprocess


#- Read in star_file
#- and sun_file
star_file = sys.argv[1]
fid = open(star_file,'r')
star_data = fid.read()

sun_file = sys.argv[2]
fid = open(sun_file,'r')
sun_data = fid.read()


#- Read value for ew_corr option
#-
ew_corr = int(sys.argv[3])

#- Read in Teff, log_g, fe_h, and vmicro
#- from the model atmosphere
#-
Teff = int(sys.argv[4])
logg = float(sys.argv[5])
feh = float(sys.argv[6])
vmicro = float(sys.argv[7])
logfile = sys.argv[8]
atm_count = int(sys.argv[9])
logfile2 = sys.argv[10]
n_Fe1 = float(sys.argv[11])
n_Fe2 = float(sys.argv[12])
star_file_mod = (sys.argv[13]) #starfile deleted if parameters are not a solution
star_log = (sys.argv[14]) #starlog deleted if parameters are not a solution
star_model = (sys.argv[15]) #starmodel deleted if parameters are not a solution


logfile3 = 'soluxns.log'
logfile4 = 'soluxns_almost.log'
fid1 = open(logfile2,'w')


#- Extract only numerical rows from star_file
#- and sun_file
#-
#- Also count the number of lines
#- measured for each element
#-
star_lines = star_data.split('\n')
sun_lines = sun_data.split('\n')

star_lambdas = [0,0]
star_eps = [0,0]
star_ews = [0,0]
star_rews = [0,0]
star_logNs = [0,0]

elem_strs = [0,0]
elem_n_lines = [n_Fe1,n_Fe2]

j = 0
for line in star_lines:

        columns = line.split()

        try:
             star_lambda = float(columns[0])
             star_ep = float(columns[1])
             star_ew = float(columns[3])
             star_rew = float(columns[4])
             star_logN = float(columns[5])

	     #elem_count = elem_count+1

        except IndexError:
	     continue
        except ValueError:
 	     continue
      


        if j < 2:

	    star_lambdas[j] = star_lambda
	    star_eps[j] = star_ep
	    star_ews[j] = star_ew
	    star_rews[j] = star_rew
	    star_logNs[j] = star_logN


        else:

	    star_lambdas.append(star_lambda)
	    star_eps.append(star_ep)
	    star_ews.append(star_ew)
	    star_rews.append(star_rew)
	    star_logNs.append(star_logN)
	
	j = j+1
      


sun_lambdas = [0,0]
sun_eps = [0,0]
sun_ews = [0,0]
sun_rews = [0,0]
sun_logNs = [0,0]

j = 0
for line in sun_lines:

        columns = line.split()

        try:
             sun_lambda = float(columns[0])
	     sun_ep = float(columns[2])
	     sun_ew = float(columns[4])
	     sun_rew = float(columns[5])
	     sun_logN = float(columns[6])

        except IndexError:
	     continue
        except ValueError:
	     continue
      

        if j < 2:

	     sun_lambdas[j] = sun_lambda
	     sun_eps[j] = sun_ep
	     sun_ews[j] = sun_ew
	     sun_rews[j] = sun_rew
	     sun_logNs[j] = sun_logN

	else:

	     sun_lambdas.append(sun_lambda)
	     sun_eps.append(sun_ep)
	     sun_ews.append(sun_ew)
	     sun_rews.append(sun_rew)
	     sun_logNs.append(sun_logN)
	
	
	j = j+1
      


if len(star_lambdas) != len(sun_lambdas):

	print '\n ERROR: star_file and sun_file are not the same length! \n'
	exit()



#- Extract lambda, ep, ew, rew, N_Fe1, N_Fe2
#- from the star_file and the sun_file
#-

star_Fe1s = [0,0]
star_Fe1_avg = 0
star_Fe1_ep_avg = 0
star_Fe1_ew_avg = 0
star_Fe1_rew_avg = 0
star_Fe2_avg = 0
star_Fe2_ep_avg = 0
star_Fe2_ew_avg = 0
i = 0
j = 0 #say what you're initiating for each of these
for elem in star_lambdas:
    

	  if j < elem_n_lines[0]:

	      star_Fe1_lambda = star_lambdas[j]
	      star_Fe1_ep = star_eps[j]
	      star_Fe1_ew = star_ews[j]
	      star_Fe1_rew = star_rews[j]
	      star_N_Fe1 = star_logNs[j]

	      sun_Fe1_lambda = sun_lambdas[j]
	      sun_Fe1_ep = sun_eps[j]
	      sun_Fe1_ew = sun_ews[j]
	      sun_Fe1_rew = sun_rews[j]
	      sun_N_Fe1 = sun_logNs[j]
              
	      #-compute sums for avgs of [Fe1/H],
	      #-Fe1_ep,Fe1_rew for correlation
	      #-coefficients
              #-
              star_Fe1 = star_N_Fe1-sun_N_Fe1

	      if i < 2:

		 star_Fe1s[i] = star_Fe1

	      else:

		 star_Fe1s.append(star_Fe1)

	      star_Fe1_avg = star_Fe1_avg + star_Fe1
	      star_Fe1_ep_avg = star_Fe1_ep_avg + star_Fe1_ep
	      star_Fe1_ew_avg = star_Fe1_ew_avg + star_Fe1_ew
	      star_Fe1_rew_avg = star_Fe1_rew_avg + star_Fe1_rew
              
	      lambda1=star_Fe1_lambda
	      lambda2=sun_Fe1_lambda
	      logN1=star_N_Fe1
	      logN2=sun_N_Fe1

              if lambda1 != lambda2:
	         print '\nERROR! The linelists for the two stars are not identical\n'
	         exit()


              if (logfile2 != 'foo.log'):
		 fid1=open(logfile2,'a')
                 #-print format for python 2.5 or lower
	         form='%7.2f  %9.2f %7.2f %7.2f %7.2f\n' 
	         fid1.write(form % (lambda1,lambda2,logN1,logN2,star_Fe1))

          else:

	      star_Fe2_lambda = star_lambdas[j]
	      star_Fe2_ep = star_eps[j]
	      star_Fe2_ew = star_ews[j]
	      star_Fe2_rew = star_rews[j]
	      star_N_Fe2 = star_logNs[j]

	      sun_Fe2_lambda = sun_lambdas[j]
	      sun_Fe2_ep = sun_eps[j]
	      sun_Fe2_ew = sun_ews[j]
	      sun_Fe2_rew = sun_rews[j]
	      sun_N_Fe2 = sun_logNs[j]

              star_Fe2 = star_N_Fe2-sun_N_Fe2
	      star_Fe2_avg = star_Fe2_avg + star_Fe2
	      star_Fe2_ep_avg = star_Fe2_ep_avg + star_Fe2_ep
	      star_Fe2_ew_avg = star_Fe2_ew_avg + star_Fe2_ew

	      lambda1=star_Fe2_lambda
	      lambda2=sun_Fe2_lambda
	      logN1=star_N_Fe2
	      logN2=sun_N_Fe2

              if (logfile2 != 'foo.log'):
	          form='%7.2f  %9.2f %7.2f %7.2f %7.2f\n' 
	          fid1.write(form % (lambda1,lambda2,logN1,logN2,star_Fe2))

        
          j=j+1
          i=i+1


star_Fe1_avg = star_Fe1_avg/elem_n_lines[0]
star_Fe2_avg = star_Fe2_avg/elem_n_lines[1]
mean_Fe1_ep = star_Fe1_ep_avg/elem_n_lines[0]
mean_Fe1_ew = star_Fe1_ew_avg/elem_n_lines[0]
mean_Fe1_rew = star_Fe1_rew_avg/elem_n_lines[0]
mean_Fe2_ep = star_Fe2_ep_avg/elem_n_lines[1]
mean_Fe2_ew = star_Fe2_ew_avg/elem_n_lines[1]

diff = abs(star_Fe1_avg-star_Fe2_avg)

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
j=0
mean_star_Fe1=star_Fe1_avg
mean_star_Fe2=star_Fe2_avg
sum_Fe1xEP=0
sum_Fe1xREW=0
sum_REWxREW=0
sum_Fe1xFe1=0
sum_EWxEW_Fe1=0
sum_EWxEP_Fe1=0
sum_EWxEW_Fe2=0
sum_EPxEP_Fe1=0
sum_EWxEP_Fe2=0
sum_EPxEP_Fe2=0
for elem in star_ews:

      if j < elem_n_lines[0]:

             sum_EPxEP_Fe1 = sum_EPxEP_Fe1 + (star_eps[j]-mean_Fe1_ep)**2
             sum_Fe1xEP = sum_Fe1xEP + (star_eps[j]-mean_Fe1_ep)*(star_Fe1s[j]-mean_star_Fe1)
             sum_REWxREW = sum_REWxREW + (star_rews[j]-mean_Fe1_rew)**2
             sum_Fe1xREW = sum_Fe1xREW + (star_rews[j]-mean_Fe1_rew)*(star_Fe1s[j]-mean_star_Fe1)
             sum_Fe1xFe1 = sum_Fe1xFe1 + (star_Fe1s[j]-mean_star_Fe1)**2
             sum_EWxEW_Fe1 = sum_EWxEW_Fe1 + (star_ews[j]-mean_Fe1_ew)**2
             sum_EWxEP_Fe1 = sum_EWxEP_Fe1 + (star_ews[j]-mean_Fe1_ew)*(star_eps[j]-mean_Fe1_ep)

      else:

             sum_EPxEP_Fe2 = sum_EPxEP_Fe2 + (star_eps[j]-mean_Fe2_ep)**2
             sum_EWxEW_Fe2 = sum_EWxEW_Fe2 + (star_ews[j]-mean_Fe2_ew)**2
             sum_EWxEP_Fe2 = sum_EWxEP_Fe2 + (star_ews[j]-mean_Fe2_ew)*(star_eps[j]-mean_Fe2_ep)

      j = j+1



r_Fe1vsEP = sum_Fe1xEP/math.sqrt(sum_Fe1xFe1*sum_EPxEP_Fe1)
r_Fe1vsREW = sum_Fe1xREW/math.sqrt(sum_Fe1xFe1*sum_REWxREW)
r_EWvsEP_Fe1 = sum_EWxEP_Fe1/math.sqrt(sum_EWxEW_Fe1*sum_EPxEP_Fe1)

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
        


	

#-
#- END OF SCRIPT

#
#   Script generates sequences of polymer molecules with wished frequencies
#   for monomers. However, there can't be two maleic acid monomers in a row.
#

import numpy as np
import shutil, os, sys
from getopt import getopt

opts, args = getopt(sys.argv[1:], 'n:',longopts=['length=','rmsd=','out_filename=','pm=','ph=','polymer='])

#############################################
################  DEFAULTS  #################
#############################################

## Number of polymer molecules, mean of their length, and rmsd of length
n = 100
meanlen = 36
rmsd = 3

## Frequencies of monomers
pm0 = 0.25
ph0 = 1 - pm0

## Name of file with sequences
out_filename = 'sequences.txt'

## Type of polymer (dibma or sma)
polymer = 'sma'

#############################################
#############  END OF DEFAULTS  #############
#############################################

# Specified options
for o, a in opts:
        if o == '-n':
            n = int(a)
        if o == '--length':
            meanlen = float(a)
        if o == '--rmsd':
            rmsd = float(a)
        if o == '--pm':
            pm0 = float(a)
        if o == '--ph':
            ph0 = float(a)
        if o == '--polymer':
            polymer = str(a)
        if o == '--out_filename':
            out_filename = str(a)
         


## Check if polymer == 'sma' or 'dibma'
while True:
    if not ((polymer == 'sma') or (polymer == 'dibma')):
        polymer = raw_input('\nIncorrect \'--polymer\' option specified. "dibma" or "sma" are available. Specify once again, please.\n')
        continue
    else:
        break

if polymer == 'sma':
    hyd_mon = 'S'
elif polymer == 'dibma':
    hyd_mon = 'D'

hyd_dict = { 'S' : 'Styrol',
             'D' : 'Diisobulylene' }
             
## Normalizing occurencies
pm0=pm0/(pm0+ph0)
ph0=1 - pm0
print('\nNormalized occurencies are '+str(pm0)+' for Maleic acid and '+str(ph0)+' for '+str(hyd_dict[hyd_mon])+'.\n')

## Frequencies for HM/H-units joining. pm0 = phm/(2phm+ph) = phm/(1+phm)
ph = (1 - 2*pm0)/(1 - pm0)
phm = pm0/(1 - pm0)

## Correction for mean length of polymer molecules caused by joining SM-unit
## while standing in the last cell of desired length, thus overrunning it by 1.
meanlen = meanlen - 0.66667*phm

#############################################################################################
with open(out_filename, 'w') as f:
    for i in xrange(int(n)):
        N = 0
        mal_count = 0
        length = int(np.random.normal(meanlen,rmsd))
        seq_list = []
        seq_list.append(np.random.choice([hyd_mon, 'M'], p=[ph0, pm0]))
        N+=1
        
        while N < length:
            if seq_list[-1] == 'M':
                seq_list.append(np.random.choice([hyd_mon, hyd_mon + 'M'], p=[ph, phm]))
                if seq_list[-1] == hyd_mon + 'M':
                    N+=2
#                    mal_count+=1
                elif seq_list[-1] == hyd_mon:
                    N+=1
            elif seq_list[-1] == hyd_mon:
                seq_list.append(np.random.choice([hyd_mon, 'M'], p=[ph0, pm0]))
                if seq_list[-1] == hyd_mon:
                    N+=1
                elif seq_list[-1] == 'M':
                    N+=1
#                    mal_count+=1
            elif seq_list[-1] == hyd_mon + 'M':
                seq_list.append(np.random.choice([hyd_mon, hyd_mon + 'M'], p=[ph, phm]))    
                if seq_list[-1] == hyd_mon + 'M':
                    N+=2
#                    mal_count+=1 
                if seq_list[-1] == hyd_mon:
                    N+=1                   
        
        for item in seq_list:
            f.write("%s" % item)
#        f.write(str(N))#+" "+str(mal_count))
        f.write("\n")


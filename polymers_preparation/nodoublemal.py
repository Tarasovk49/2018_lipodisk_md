#
#   Scripts generates sequences of polymer molecules with wished frequencies
#   for monomers, in which there can't be two maleic acid monomers in a row.
#

import numpy as np
import shutil, os, sys

## Number of polymer molecules, mean of their length, and rmsd of length
n = 100
meanlen = 36
rmsd = 3
## Wished frequencies for monomers
pm0 = 0.25
ps0 = 1 - pm0
## Frequencies for SM/S-units joining. pm0 = psm/(2psm+ps) = psm/(1+psm)
ps = (1 - 2*pm0)/(1 - pm0)
psm = pm0/(1 - pm0)
## Correction for mean length of polymer molecules caused by joining SM-unit
## while standing in the last cell of desired length, thus overrunning it by 1.
meanlen = meanlen - 0.66667*psm
with open('seq.txt', 'w') as f:
    for i in xrange(int(n)):
        N = 0
        mal_count = 0
        length = int(np.random.normal(meanlen,rmsd))
        seq_list = []
        seq_list.append(np.random.choice(["S","M"], p=[ps0, pm0]))
        N+=1
        
        while N < length:
            if seq_list[-1] == "M":
                seq_list.append(np.random.choice(["S","SM"], p=[ps, psm]))
                if seq_list[-1] == "SM":
                    N+=2
#                    mal_count+=1
                elif seq_list[-1] == "S":
                    N+=1
            elif seq_list[-1] == "S":
                seq_list.append(np.random.choice(["S","M"], p=[ps0, pm0]))
                if seq_list[-1] == "S":
                    N+=1
                elif seq_list[-1] == "M":
                    N+=1
#                    mal_count+=1
            elif seq_list[-1] == "SM":
                seq_list.append(np.random.choice(["S","SM"], p=[ps, psm]))    
                if seq_list[-1] == "SM":
                    N+=2
#                    mal_count+=1 
                if seq_list[-1] == "S":
                    N+=1                   
        
        for item in seq_list:
            f.write("%s" % item)
#        f.write(str(N))#+" "+str(mal_count))
        f.write("\n")


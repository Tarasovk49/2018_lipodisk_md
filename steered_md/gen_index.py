# Generates .ndx file with groups needed for pulling: groups POL$N including backbone atoms
# of polymers and group MEMBRANE including phosphorus atoms of membrane. These groups then
# used during simulation to compute centers of mass of polymers and membrane.
# Forces are applied on distances between $MEMBRANE and POL$N groups to pull polymers towards membrane.

import sys
from getopt import getopt

opts, args = getopt(sys.argv[1:], 'i:o:')

in_filename = 'lipodisk.pdb'
out_filename = 'index_more.ndx'

for o, a in opts:
    if o == '-i':
        in_filename = a
    if o == '-o':
        out_filename = a

def make_ndx(in_filename, out_filename):
    """
    Creates gromacs index file with groups: POL@X - backbone carbons of polymer molecule X
    and MEMBRANE - phosphorus atoms of lipids.
    """
    with open(in_filename) as pdb:
        with open(out_filename,'w') as ndx:
            resid=1000
            pol_molecule_num = 1 # Counts polymer molecules
            line_value_num = 0 # Counts written numbers of atoms in strings (numbers are written by 15 in one string)
            for line in pdb.readlines():
                if (line[:4]=='ATOM') and ((line[17:21]=='POPC') or (line[17:21]=='DPPC')\
                 or (line[17:21]=='DMPC') or (line[17:21]=='DOPC') or (line[17:21]=='CHOL')) \
                                      and ('P' in line[12:16]):
                    if int(line[23:26] )<resid:
                        ndx.write('\n' + '[ MEMBRANE ]' + '\n')
                        resid = int(line[23:26])
                        line_value_num = 0
                    ndx.write(' '+line[6:12])
                    line_value_num += 1
                    if line_value_num == 15:
                        ndx.write('\n')
                        line_value_num = 0
                if (line[:4]=='ATOM') and ((line[17:20]=='ST1') or (line[17:20]=='ST2') or (line[17:20]=='MAL')\
                     or (line[17:20]=='MAR') or (line[17:20]=='MAD') or (line[17:20]=='MA2') or (line[17:20]=='DB1')\
                     or (line[17:20]=='DB2')) and (('CA' in line[12:16]) or ('CB' in line[12:16])):
                    if line_value_num == 15:
                        ndx.write('\n')
                        line_value_num = 0
                    if int(line[23:26] )<resid:
                        ndx.write('\n' + '[ POL' + str(pol_molecule_num) + ' ]' + '\n')
                        ndx.write(' '+line[6:12])
                        resid = int(line[23:26])
                        pol_molecule_num +=1
                        line_value_num = 0
                    elif int(line[23:26])==resid:
                        ndx.write(' '+line[6:12])
                    else:
                        resid+=1
                        ndx.write(' '+line[6:12])
                    line_value_num += 1
#    print pol_molecule_num
    
make_ndx(in_filename, out_filename)

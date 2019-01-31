import numpy as np
import MDAnalysis as mda
from sys import argv
from getopt import getopt
import os

opts, args = getopt(argv[1:], 'm:o:')

membrane = 'popc.pdb'
outfile = 'lipodisk.pdb'

for o, a in opts:
    if o == '-m':
        membrane = a
    if o == '-o':
        outfile = a

# Load membrane
u_mem = mda.Universe(membrane)
ag_mem = u_mem.select_atoms("resname POPC DOPC CHOL DPPC")
# Place membrane in (0,0,0)
ag_mem.translate(-ag_mem.center_of_geometry())
# Compute radius of membrane
xmax,ymax,zmax=np.amax(ag_mem.positions[:],0)
xmin,ymin,zmin=np.amin(ag_mem.positions[:],0)
R = np.sqrt((xmax-xmin)**2 + (ymax-ymin)**2)/2


# Number of polymer layers
layers = 3
# Decrease 2 coefficients below to place polymers with higher density
# r_coeff = 1.0 is the least value when it is garanteed that none of the polymers overlap
r_coeff = 0.7
# alpha_coeff = 2.0 is the least value when it is garanteed that none of the polymers overlap
angle_coeff = 0.9


# Load polymers, place them in (0,0,0), compute the maximal radius
pol_list = os.listdir('initial_structures')
polymers={}
rmax = 0
for pol_num, pol_name in enumerate(pol_list):
    u = mda.Universe('initial_structures/%s'%pol_name)
    polymers[pol_num] = u.select_atoms("resname DB1 DB2 MAR MAL MA2 MAD ST1 ST2")
    polymers[pol_num].translate(-polymers[pol_num].center_of_geometry())
    PAX = np.matrix(polymers[pol_num].principal_axes())
    PAXI = PAX.I
    polymers[pol_num].rotate(PAX)
    print(polymers[pol_num].principal_axes())
    xmax, ymax, zmax = np.amax(ag_mem.positions[:],0)
    xmin, ymin, zmin = np.amin(ag_mem.positions[:],0)
    r = np.sqrt((xmax-xmin)**2 + (ymax-ymin)**2 + (zmax-zmin)**2)/2
    if r > rmax:
        rmax = r

# Place polymers around membrane
lastvalue = 0
for b in range(layers):
#    b = b + 1
    angle = angle_coeff*np.arcsin(rmax/(rmax*(2*b+1)*r_coeff + R))
    k = int(6.28318530718/angle)
    angle = 6.28318530718/k
    for a in range(k):
        pol_num = lastvalue + a
        x = np.cos(angle*pol_num)*(rmax*(2*b+1)*r_coeff + R)
        y = np.sin(angle*pol_num)*(rmax*(2*b+1)*r_coeff + R)
        polymers[pol_num].translate([x, y, 0.0])
        u_mem = mda.core.universe.Merge(u_mem.atoms, polymers[pol_num].atoms)
    lastvalue = pol_num + 1

u_mem.atoms.write(outfile)

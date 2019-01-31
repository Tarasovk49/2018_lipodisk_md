import numpy as np
import MDAnalysis as mda
from sys import argv
from getopt import getopt
import os

opts, args = getopt(argv[1:], 'm:p:o:n:')

# Defaults
membrane = 'popc.pdb'
pol_dir = 'polymer_molecules'
outfile = 'lipodisk.pdb'
layers = 2

for o, a in opts:
    if o == '-m':
        membrane = a
    if o == '-p':
        pol_dir = a
    if o == '-o':
        outfile = a
    if o == '-n':
        layers = a
        
# Load membrane
u_mem = mda.Universe(membrane)
ag_mem = u_mem.select_atoms("resname POPC DOPC CHOL DPPC")
# Place membrane in (0,0,0)
ag_mem.translate(-ag_mem.center_of_geometry())
# Compute radius of membrane
xmax,ymax,zmax=np.amax(ag_mem.positions[:],0)
xmin,ymin,zmin=np.amin(ag_mem.positions[:],0)
R = np.sqrt((xmax-xmin)**2 + (ymax-ymin)**2)/2
print('R of membrane is '+R)

# Decrease 2 coefficients below to place polymers with higher density
# r_coeff = 1.0 is the least value when it is garanteed that none of the polymers overlap
r_coeff = 0.7
# angle_coeff = 2.0 is the least value when it is garanteed that none of the polymers overlap
angle_coeff = 0.9


# Load polymers, place them in (0,0,0), compute the maximal radius
pol_list = os.listdir(pol_dir)
polymers={}
rmax = 0
for pol_num, pol_name in enumerate(pol_list):
    u = mda.Universe(pol_dir + '/%s'%pol_name)
    polymers[pol_num] = u.select_atoms("resname DB1 DB2 MAR MAL MA2 MAD ST1 ST2")
    polymers[pol_num].translate(-polymers[pol_num].center_of_geometry())
    PAX = np.matrix(polymers[pol_num].principal_axes())
    PAXI = PAX.I
    polymers[pol_num].rotate(PAX)
    xmax, ymax, zmax = np.amax(ag_mem.positions[:],0)
    xmin, ymin, zmin = np.amin(ag_mem.positions[:],0)
    r = np.sqrt((xmax-xmin)**2 + (ymax-ymin)**2)/2
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

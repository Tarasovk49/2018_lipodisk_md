#
#             Calculation of lipid order parameters.
#
# DESCRIPTION
# Script follows next steps:
#
#   1. Fits plane to selected set of dots (C27 and C37 atoms of lipids) - selection_planefit
#   2. Computes order parameters for each subset of lipid atoms across normal to computed plane.
#      Subsets include (13 lipid fatty acid chains carbon atoms - C2, C3, .., C14) * (3 subsets of lipids -
#      all, inner, outer) = 39 subsets
#
# Output includes:
#    - plot Order parameter vs atom number in lipid fatty acid chain


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
from tqdm import notebook
import scipy.optimize
from mpl_toolkits.mplot3d import Axes3D
import os, sys
from getopt import getopt

opts, args = getopt(sys.argv[1:], 's:e:i:', longopts=['step=','write_freq=','structure=','trajectory=','planefit_sel='])


#############################################
################  DEFAULTS  #################
#############################################

step = 2 # Step of integration in fs
write_freq = 5000 # Number of steps that elapse between writing coordinates to output trajectory file
interval = 1 # Interval to compute new normal to membrane plane (0.5 ns)
s = 0 # in ps
e = 4000 # in ps

membrane = 'lipodisk_npt.gro'
traj = 'lipodisk_npt_2_whole_cluster_nojump_mol.xtc'

selection_planefit = 'not resid 111 106 and resname DMPC and name C27 C37'

#############################################
#############  END OF DEFAULTS  #############
#############################################


# Specified options
for o, a in opts:
        if o == '-s':
            start_frame = int(a)
        if o == '-e':
            end_frame = int(a)
        if o == '-i':
            interval = int(a)
        if o == '--step':
            step = int(a)
        if o == '--write_freq':
            write_freq = int(a)
        if o == '--structure':
            membrane = str(a)
        if o == '--trajectory':
            traj = str(a)
        if o == '--planefit_sel':
            selection_planefit = str(a)

# Transform picoseconds to frames
start_frame = int(s/step)
end_frame = int(e/step)
interval = int(500000/(step*write_freq))

universe = mda.Universe(membrane, traj)
plane_group = universe.select_atoms(selection_planefit)



def fitPlaneLTSQ(XYZ):
    '''
    Fits plane to set of points
    Source:
    https://gist.github.com/RustingSword/e22a11e1d391f2ab1f2c
    '''
    (rows, cols) = XYZ.shape
    G = np.ones((rows, 3))
    G[:, 0] = XYZ[:, 0]  #X
    G[:, 1] = XYZ[:, 1]  #Y
    Z = XYZ[:, 2]
    (a, b, c),resid,rank,s = np.linalg.lstsq(G, Z)
    normal = (a, b, -1)
    nn = np.linalg.norm(normal)
    normal = normal / nn
    return (c, normal)

def calculate_S(coordinates, normal):
    '''
    Calculates order parameter S for set of coordinates across normal
    '''
    cd = np.concatenate((coordinates[1::3]-coordinates[0::3], coordinates[2::3]-coordinates[0::3]), axis=0)
    cd_r = np.sqrt(np.sum(np.power(cd,2), axis=1))
    matrix_multiply = np.matmul(cd,normal)
    cos_theta = matrix_multiply/cd_r
    S_cd = -0.5*(3.*np.square(cos_theta)-1)
    order_param = np.average(S_cd)
    return order_param


LOP_group_all = {}
LOP_group_inner = {}
LOP_group_outer = {}
for atom_number in ['2','3','4','5','6','7','8','9','10','11','12','13','14']:
    selection_all = 'not resid 111 106 and resname DMPC and (name C2{0} or name H{0}R or name H{0}S or name C3{0} or name H{0}X or name H{0}Y)'.format(atom_number)
    selection_inner = 'not resid 111 106 and resname DMPC and (name C2{0} or name H{0}R or name H{0}S or name C3{0} or name H{0}X or name H{0}Y) and not same resid as (around 10 resname MAL MAR ST1 ST2 DB1 DB2 MAD MA2)'.format(atom_number)
    selection_outer = 'not resid 111 106 and resname DMPC and (name C2{0} or name H{0}R or name H{0}S or name C3{0} or name H{0}X or name H{0}Y) and same resid as (around 10 resname MAL MAR ST1 ST2 DB1 DB2 MAD MA2)'.format(atom_number)
    LOP_group_all[atom_number] = universe.select_atoms(selection_all)
    LOP_group_inner[atom_number] = universe.select_atoms(selection_inner)
    LOP_group_outer[atom_number] = universe.select_atoms(selection_outer)
    
    
data_all = pd.DataFrame(columns = ['C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14'])
data_inner = pd.DataFrame(columns = ['C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14'])
data_outer = pd.DataFrame(columns = ['C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14'])

for ts in notebook.tqdm(universe.trajectory):
    if (ts.frame >= start_frame) and (ts.frame <= end_frame):
        if ts.frame%interval == 0:
            positions = plane_group.atoms.positions[:]
            c, normal = fitPlaneLTSQ(positions)

        OP_list_all = []
        OP_list_inner = []
        OP_list_outer = []
        for atom_number in ['2','3','4','5','6','7','8','9','10','11','12','13','14']:
            # ALL
            coordinates_all = LOP_group_all[atom_number].atoms.positions[:]
            order_param_all = calculate_S(coordinates_all,normal)
            OP_list_all.append(order_param_all)
            # INNER
            coordinates_inner = LOP_group_inner[atom_number].atoms.positions[:]
            order_param_inner = calculate_S(coordinates_inner,normal)
            OP_list_inner.append(order_param_inner)
            # OUTER
            coordinates_outer = LOP_group_outer[atom_number].atoms.positions[:]
            order_param_outer = calculate_S(coordinates_outer,normal)
            OP_list_outer.append(order_param_outer)
        # ALL    
        a_series_all = pd.Series(OP_list_all, index = data_all.columns)
        data_all = data_all.append(a_series_all, ignore_index=True)
        # INNER
        a_series_inner = pd.Series(OP_list_inner, index = data_inner.columns)
        data_inner = data_inner.append(a_series_inner, ignore_index=True)
        # OUTER
        a_series_outer = pd.Series(OP_list_outer, index = data_outer.columns)
        data_outer = data_outer.append(a_series_outer, ignore_index=True)
    elif ts.frame > end_frame:
        break
        
plt.figure(figsize=(10, 7), dpi=80)
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)\
data_all.mean().rename('All').plot(legend='True')
data_inner.mean().rename('Inner').plot(legend='True')
data_outer.mean().rename('Outer').plot(legend='True')
plt.ylabel(u\"Order parameter, S\", fontsize=14),
plt.xlabel(u\"Fatty acid chain carbon\", fontsize=14)
plt.grid()
#plt.show()
plt.savefig('Lipid_OP_dynamics.png')

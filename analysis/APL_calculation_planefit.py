#
#             Calculation of area per lipid.
#
# DESCRIPTION
# Script follows next steps:
#
#   1. Fits plane to selected set of dots (C27 atoms of lipids) - selection1
#   2. Finds projections onto that plane of selection1 and selection2
#      representing another set of dots (Protein CA atoms and polymer CA CB atoms)
#   3. Gets rid of protein and polymer dots hanging over membrane by deleting from selection2 atoms
#      which distances between atom positions and their projections are greater than 20 angstroms
#   4. Rotates projections in such a way that they now lie in XY-plane
#   5. Builds Delaunay triangulation for these points
#   6. Computes areas of triangles
#   7. Computes SUM of areas of triangles in such a way:
#     - triangles formed by only lipid atoms contribute their whole area to SUM
#     - triangles formed by lipid and polymer atoms contribute half of their area to SUM
#     - triangles formed by only polymer atoms don't contribute to SUM
#   8. SUM is then divided by number of lipids and multiplied by 2 cause both monolayers
#      were accounted for when computing projections
#  
#
#
# Output includes:
#     - CSV file with (time_step, APL)
#     - plot APL vs time_step
#     - 3D plot of prefered sets of points obtained on different stages of script
#       (to check if script works properly)
#     - Delaunay triangulation plot

from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
import math
import MDAnalysis as mda
import mpld3
from mpl_toolkits.mplot3d import Axes3D
import os, sys
from getopt import getopt


opts, args = getopt(sys.argv[1:], 's:e:i:',longopts=['structure=','trajectory=','planefit_sel=','dots_sel='])

#############################################
################  DEFAULTS  #################
#############################################

## Start and end frames to calculate APL, interval to recalculate plane
start_frame = 0
end_frame = 1000000
interval = 1


#membrane = 'SMALP_npt.gro'
membrane = 'SMALP_ions.gro'
#traj = 'SMALP_npt.xtc'
traj = 'test.xtc'



selection1 = 'resname DMPC and name C27'
selection2 = '(resname ST1 ST2 MAL MAR DB1 DB2 MAD MA2 and name CA CB) or (protein and name CA)'

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
        if o == '--structure':
            membrane = str(a)
        if o == '--trajectory':
            traj = str(a)
        if o == '--planefit_sel':
            selection1 = str(a)
        if o == '--dots_sel':
            selection2 = str(a)


universe = mda.Universe(membrane, traj)



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

def classify( x ):
    '''
    Classifies the triangle of Delaunay triangulation
    as triangle formed by lipids (a=1), by lipids and polymers (a=2)
    or by polymers alone (a=3)
    '''
    a = 1
    lipidexist = False
    polexist = False
    for i in x:
        if i < border:
            lipidexist = True
        if i >= border:
            polexist = True
    if lipidexist and polexist:
        a = 2
    elif not lipidexist and polexist:
        a = 3
    return a

def rotation_matrix(axis, theta):
    """
    Returns the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    Source:
    https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
                     
                     
def get_projected_points( normal, points ):
    '''
    Finds projections of points onto a plane given by normal
    with the point (0,0,0). Then rotates those projections into xy plane.
    Source:
    https://math.stackexchange.com/questions/100761/how-do-i-find-the-projection-of-a-point-onto-a-plane
    https://math.stackexchange.com/questions/1167717/transform-a-plane-to-the-xy-plane
    '''
    # Calculate array of displaces for each point to get projections onto plane with @normal and point (0,0,0)
    t = - np.einsum('j,ij', normal, points)
    # Apply those displaces to each point to get projections of points onto plane with @normal and point (0,0,0)
    projected_points = points + np.outer(t,normal)
    # Derive angle between plane with @normal and xy-plane with (0,0,1)
    cos_theta = np.dot(normal,[0,0,1])
    theta = np.arccos(cos_theta)
    # Derive axis of rotation
    axis = np.cross(normal,[0,0,1])
    # Rotate points' projections onto xy-plane
    final_points_3D = np.matmul(rotation_matrix(axis, theta), projected_points.T).T
    # Reduce 3D-array to 2D-array of points' coordinates by deleting z-coordinate
    final_points_2D = np.delete(final_points_3D,obj=2,axis=1)
    # You can return all sets of points to plot them and check if script works properly
    return final_points_2D
    #return final_points_2D, projected_points, final_points_3D
    
def get_projected_points_1(c, normal, points ):
    '''
    Finds projections of points onto a plane given by normal
    with the point (0,0,c).
    Source:
    https://math.stackexchange.com/questions/100761/how-do-i-find-the-projection-of-a-point-onto-a-plane
    '''
    # Calculate array of displaces for each point to get projections onto plane with @normal and point (0,0,c)
    subtract = points - np.array([0,0,c])
    t = np.einsum('j,ij', normal, subtract)
    # Apply those displaces to each point to get projections of points onto plane with @normal and point (0,0,c)
    projected_points = points + np.outer(t,normal)
    return projected_points

def get_projected_points_2( normal, points ):
    '''
    Rotates projections into xy plane.
    Source:
    https://math.stackexchange.com/questions/1167717/transform-a-plane-to-the-xy-plane
    '''
    # Derive angle between plane with @normal and xy-plane with (0,0,1)
    cos_theta = np.dot(normal,[0,0,1])
    theta = np.arccos(cos_theta)
    # Derive axis of rotation
    axis = np.cross(normal,[0,0,1])
    # Rotate points' projections onto xy-plane
    final_points_3D = np.matmul(rotation_matrix(axis, theta), points.T).T
    # Reduce 3D-array to 2D-array of points' coordinates by deleting z-coordinate
    final_points_2D = np.delete(final_points_3D,obj=2,axis=1)
    # You can return all sets of points to plot them and check if script works properly
    return final_points_2D
    #return final_points_2D, projected_points, final_points_3D
    

lipid_group = universe.select_atoms(selection1)
polymer_group = universe.select_atoms(selection2)



ts_list = []
apl_list = []
for ts in universe.trajectory:
    if (ts.frame >= start_frame) and (ts.frame <= end_frame):
        if ts.frame%interval == 0:
            # Fit plane to lipid atoms(selection1)
            positions = lipid_group.atoms.positions[:]
            c, normal = fitPlaneLTSQ(positions)
        # Get atom coordinates
        lipid_atoms = lipid_group.atoms.positions[:]
        polymer_atoms = polymer_group.atoms.positions[:]
        # Get 2D atom projections onto plane fitted to lipid headgroups
        points_lipid = get_projected_points(normal, lipid_atoms)
        # There two options to do that:
        #       1. Count all polymer atoms including those which are hanging over lipid membrane
        #points_pol = get_projected_points(normal, polymer_atoms)
        
        #       2. Doesn't count polymer atoms which are hanging over lipid membrane
        #         (distance between their positions and their projections are greater than 20 angstroms)
        points_pol_1 = get_projected_points_1(c, normal, polymer_atoms)
        pol_delta = polymer_atoms-points_pol_1
        dist = np.sqrt(np.einsum('ij,ij->i', pol_delta, pol_delta))
        points_pol_1 = points_pol_1[dist<20]
        points_pol_2 = get_projected_points_2(normal, points_pol_1)
        points_pol = points_pol_2
        
        #points_lipid, projected_points, rotated_points = get_projected_points(normal, lipid_atoms)
        #points_pol, ghgh, fhf = get_projected_points(normal, polymer_atoms)
        points = np.concatenate((points_lipid, points_pol), axis=0)
        # Plot Delaunay tesselation
        dela = Delaunay(points)
        # border = the first index in final_array which stands for polymer record = number of lipids
        border = points_lipid.shape[0]
        simplices_array = np.column_stack((dela.simplices[:], np.apply_along_axis(classify, axis=1, arr=dela.simplices[:])))
        SUM = 0
        # Calculate Area
        for triangle in simplices_array:
            if triangle[3]==3:
                pass
            elif triangle[3]==2:
                pass
                Ax, Ay = points[triangle[0],0], points[triangle[0],1]
                Bx, By = points[triangle[1],0], points[triangle[1],1]
                Cx, Cy = points[triangle[2],0], points[triangle[2],1]
                SUM+=abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/4
            elif triangle[3]==1:
                Ax, Ay = points[triangle[0],0], points[triangle[0],1]
                Bx, By = points[triangle[1],0], points[triangle[1],1]
                Cx, Cy = points[triangle[2],0], points[triangle[2],1]
                SUM+=abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2
        APL = SUM/border
        apl_list.append(APL)
        ts_list.append(ts.frame)
    else:
        break
        
# APL must be multiplied by 2 because lipids from both monolayers were accounted for
# Also if used not single atom of each lipid molecule
# APL is to be multiplied by number of atoms used to get correct APL values
new_apl_list = np.multiply(apl_list,2)


# Plot 3D of different sets of points 
fig = plt.figure(num=None, figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')
ax = fig.gca(projection='3d')

#selection3 = 'resname DMPC and name C28'
#test_group = universe.select_atoms(selection3)
#test_atoms = test_group.atoms.positions[:]
#ax.scatter(test_atoms[:, 0], test_atoms[:, 1], test_atoms[:, 2], color='g')

ax.scatter(lipid_atoms[:, 0], lipid_atoms[:, 1], lipid_atoms[:, 2], color='r')
#ax.scatter(projected_points[:, 0], projected_points[:, 1], projected_points[:, 2], color='b')
#ax.scatter(rotated_points[:, 0], rotated_points[:, 1], rotated_points[:, 2], color='g')
#ax.scatter(points_lipid[:, 0], points_lipid[:, 1], np.zeros(points_lipid.shape[0]), color='black')

ax.scatter(polymer_atoms[:, 0], polymer_atoms[:, 1], polymer_atoms[:, 2], color='orange')

# plot fitted plane
maxx = np.max(lipid_atoms[:,0]+100)
maxy = np.max(lipid_atoms[:,1]+100)
minx = np.min(lipid_atoms[:,0]-100)
miny = np.min(lipid_atoms[:,1]-100)

point = np.array([0.0, 0.0, c])
d = -point.dot(normal)

# compute needed points for plane plotting
xx, yy = np.meshgrid([minx, maxx], [miny, maxy])
z = (-normal[0]*xx - normal[1]*yy - d)*1. / normal[2]

# plot plane
ax.plot_surface(xx, yy, z, alpha=0.2)

ax.set_zlim3d(-110, 250)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
#fig.savefig('3D_sets_of_points.png')


# Save (time, APL) as CSV file
np.savetxt('APL_vs_time.csv', np.c_[ts_list, new_apl_list], fmt=('%.d','%.4f'), delimiter='   ', header='APL vs frame. 13 nm DMPC-SMALP with 0.0 charge of maleic acid monomer')

# Plot APL vs time
plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
plt.plot(ts_list, new_apl_list)
#plt.savefig('APL_vs_time.png')

# Plot Delaunay triangulation
plt.figure(num=None, figsize=(13, 13), dpi=80, facecolor='w', edgecolor='k')
plt.triplot(points[:,0], points[:,1], dela.simplices)
plt.plot(points_lipid[:,0], points_lipid[:,1], 'o', color='r')
plt.plot(points_pol[:,0], points_pol[:,1], 'o', color='g')
plt.savefig('Delaunay.png')                    

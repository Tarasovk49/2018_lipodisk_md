# Lab journal


### 02.08.19
#### Completely rewritten [`polymer_generation.py`](polymers_preparation/polymer_generation.py).
It now has friendly inteface and also there are three ways to build polymer molecules: 
1. Generate random polymer molecules.
2. Generate polymer molecules from txt file.
3. Generate polymer molecules from manually submitted sequences.

### 01.08.19
#### Analysis of trajectory 1h2s in DMPC layer of ~11\*11 nm^2. Whole system size is 20\*20\*12 nm^3.
First we need to generate index file with groups needed for analysis:
```
gmx_2018 make_ndx -f lipodisk_npt.tpr<<EOF
a C211 | a C212 | a C213
name 28 C11-C12-C13-myristoil
r MAL | r MAR | r MA2 | r MAD
a CA | a CB
29 & 30
name 31 MA_backbone
r ST1 | r ST2
32 & 30
name 33 ST_backbone
33 | 31
name 34 both_backbone
!27
name 35 except_wat_and_ions
q
EOF
```

C11-C12-C13-myristoil is a group for ACF calculation

MA_backbone, ST_backbone, both_backbone are groups for gyration radii calculation

except_wat_and_ions is a group for SASA calculation

**Plots:**
1. Lipid density versus thickness.
```
gmx_2018 density -f lipodisk_npt.xtc -s lipodisk_npt.tpr -n index.ndx -o lipid_density.xvg -dt 10000 -center -d Z#<<!
DMPC
DMPC
!
python density_profile.py
```
<p align="center">
  <img width="500" height="350" src="images/lipid_density.png">
</p>

2. SASA versus time.
```
gmx_2018 sasa -f lipodisk_npt.xtc -n index.ndx -s lipodisk_npt.tpr -o sasa.xvg<<!
except_wat_and_ions
!
python sasa.py
```
<p align="center">
  <img width="500" height="350" src="images/sasa.png">
</p>

3. Gyration radii of Maleic acid monomers, Styrol monomers and both versus time.
```
gmx_2018 gyrate -f lipodisk_npt.xtc -s lipodisk_npt.tpr -n index.ndx -o MA_gyr.xvg<<!
MA_backbone
!
gmx_2018 gyrate -f lipodisk_npt.xtc -s lipodisk_npt.tpr -n index.ndx -o ST_gyr.xvg<<!
ST_backbone
!
gmx_2018 gyrate -f lipodisk_npt.xtc -s lipodisk_npt.tpr -n index.ndx -o both_gyr.xvg<<!
both_backbone
!
python gyration_radii.py
```
<p align="center">
  <img width="500" height="350" src="images/gyr_radii.png">
</p>

4. Rotational ACF of myristoil atoms C11 C12 C13 (as template POPC where C15 C16 C17 were chosen in previous calculations).
```
gmx_2018 rotacf -f lipodisk_npt.xtc -s lipodisk_npt.tpr -n index.ndx -o lipodisk_rotacf.xvg -P 2<<!
C11-C12-C13-myristoil
!
python acf_fit.py
```
<p align="center">
  <img width="1000" height="300" src="images/lipodisk_rotacf.png">
</p>


#### Uploaded 2 complete algorithms
1. Charge dependence simulation.
2. Sensory rhodopsin 1h2s preparation and simulation from the very beginning.

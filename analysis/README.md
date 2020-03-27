# Description

This directory provides analysis scripts of Gromacs trajectories for lipodisks with example outputs.

### Calculate area per lipid
Calculation of area per lipid is provided by two different scripts.
1. The first one works much faster and calculates area per lipid for projections of lipids onto a fitted to those lipids plane - [*APL_calculation_planefit.py*](APL_calculation_planefit.py).

`APL_calculation_planefit.py -s (0) --e (1000000) --i (1) --structure (SMALP_ions.gro) --trajectory (SMALP_npt.xtc) --planefit_sel (resname DMPC and name C27) --dots_sel ((resname ST1 ST2 MAL MAR DB1 DB2 MAD MA2 and name CA CB) or (protein and name CA))`

- `-s 0` start frame

- `-e 1000000` end frame

- `-i 1` interval to update fitted plane

- `--structure SMALP_ions.gro` PDB or GRO file

- `--trajectory SMALP_npt.xtc` XTC or TRR file

- `--planefit_sel 'resname DMPC and name C27'` selection of atoms to fit plane

- `--dots_sel '(resname ST1 ST2 MAL MAR DB1 DB2 MAD MA2 and name CA CB) or (protein and name CA)'` additional to `--planefit_sel` selection of atoms to calculate Delaunay triangulation

2. The other one finds projections of lipids onto a fitted to those lipids parabola. That is useful if your membrane bended during the simulation. It is quite normal for membranes with area greater than ~200 nm^2 - [*APL_calculation_parabolafit.py*](APL_calculation_parabolafit.py).

It also contains original Jupyter notebooks fro those two scripts - [*APL_calculation_planefit.ipynb*](APL_calculation_planefit.ipynb), [...](...).
### Calculate order parameters of lipids

[*LOP_all_atoms.py*](LOP_all_atoms.py) calculates lipid order parameters for each type of lipid fatty acid chain carbon atom - C2, C3, C4, C5, ..., C14. Order parameters are also calculated for three different subsets of lipids - for all lipids, for central (inner) lipids that located at distance more than 2 nm from polymer molecules and peripheral (outer) lipids of lipodisk that located at distance less than 2 nm from polymer molecules.

[Original Jupyter notebook](LOP_all_atoms.ipynb) is also available.

### Calculate rotational correlation time of lipids

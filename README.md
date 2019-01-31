# Polybuild

**Polybuild** contains scripts and files to perform GROMACS simulations of lipodisks. Lipodisk is extracting and delivering system for membrane proteins. Lipodisks are comprised of a part lipid membrane, stabilized on ridges by amphiphilic polymers. The polymers usually used for producing lipodisks are SMA (Styrene and Maleic Acid monomers) and DIBMA (Diisobutylene and Maleic Acid monomers). **Polybuild** contains:

1. Script for building desired polymer molecules consisting of 3 types of monomers: Maleic acid, Diisobutylene, Styrene.
2. OPLS-AA forcefield for generating topology with GROMACS command pdb2gmx for systems with lipids and polymers.
3. Script for producing starting configuration of lipodisk for GROMACS simulations.
4. Scripts to prepare all necessary files for steered GROMACS simulation to move polymer molecules towards membrane to form lipodisk.

**Dependencies:**
1. python 3.*
2. [Pysimm](https://github.com/Tarasovk49/pysimm) - branch contains changes providing ability to read atom names from *.mol2* files and write them and also resids and resnames to *.pdb* files.
3. [MDAnalysis](https://github.com/MDAnalysis/mdanalysis).
4. [Gromacs 2018.1](http://www.gromacs.org/)

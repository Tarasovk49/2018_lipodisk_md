# Lipodisk_md

### Overview
**Lipodisk_md** contains scripts and files to perform GROMACS simulations of lipodisks. Lipodisk is membrane system for extracting membrane proteins. Lipodisk comprises a part of lipid membrane, stabilized on ridges by amphiphilic polymers. The polymers usually used for producing lipodisks are SMA (Styrene and Maleic Acid copolymer) and DIBMA (Diisobutylene and Maleic Acid copolymer). Lipodisks are therefore called SMALPs and DIBMALPs. **Lipodisk_md** allows:

1. Build polymer molecules using [Pysimm](https://github.com/Tarasovk49/pysimm) consisting of 3 types of monomers: Maleic acid, Diisobutylene, Styrene.
![lol](images/monomers.jpg)

<p align="center">
  <img width="600" height="300" src="images/monomers.jpg">
</p>
2. Generate topology files with command *pdb2gmx* for systems, containing 4 types of lipids (POPC, DOPC, DPPC, CHOL) and polymer molecules. Topologies for lipids were taken from [1](https://www.ncbi.nlm.nih.gov/pubmed/26568975), [2](https://www.ncbi.nlm.nih.gov/pubmed/24745688), [3](https://www.ncbi.nlm.nih.gov/pubmed/26187855).
3. Produce starting configuration of lipodisk via [MDAnalysis](https://github.com/MDAnalysis/mdanalysis) for [Gromacs](http://www.gromacs.org/) simulations.
4. Prepare all necessary files for steered molecular dynamics to move polymer molecules towards membrane.

### Dependencies
1. python 3.*
2. [Pysimm](https://github.com/Tarasovk49/pysimm) - branch contains changes providing ability to read atom names from *.mol2* files and write them and also resids and resnames to *.pdb* files. LAMMPS installation is also needed (see Pysimm [installation guide](https://github.com/Tarasovk49/pysimm#complete-installation-pysimm-and-lammps)).
3. [MDAnalysis](https://github.com/MDAnalysis/mdanalysis).
4. [Gromacs 2018.1.](http://manual.gromacs.org/documentation/)

### NOTE
For all Gromacs simulations use the forcefield provided (*oplsaa_lipids_polymers.ff* and *residuetypes.dat*). It includes all necessary files for polymers and for 4 types of lipids.

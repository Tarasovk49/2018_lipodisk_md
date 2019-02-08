# Description

*polymer_generation.py* builds polymer molecules consisting of diisobutylene, maleic acid, styrene monomers. 8 types of monomers are built with structure provided by *topology_generalized/.mol*:
- 4 types of maleic acid monomers = (2 different states of protonation) x (2 different carboxyls of maleic acid),
- 2 types of diisobutylene monomers with radical group oriented closer to beginning or closer to end of polymer,
- 2 types of styrene monomers with benzene group oriented closer to beginning or closer to end of polymer.

The syntax is (default values are displayed):

`polymer_generation.py -n 100 -l 36 -r 3 -p 3 --ps 0 --pm 0.25 --pd 0.75 -o polymer`

Options are:

**-n** Number of polymer molecules to generate

**-meanlen** Mean length of polymer molecules

**-rmsd** Root mean square deviation of length of polymer molecules

**-protstate** Protonation state is a number that defines the mean charge of maleic acid monomer. Correspondence of protonation state, pH value and mean charge of MA monomer is given by table:

| Protonation state   | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| pH  | 5 | 6 | 7 | 8 | 9 | 10 | >10 |
| charge | -0.3 | -0.5 | -1.0 | -1.2 | -1.7 | -1.9 | -2.0 |

**-pm** Occurrence of maleic acid monomer

**-pd** Occurrence of diisobutylene monomer

**-ps** Occurrence of styrene monomer

**-out_filename** Name of generated *.pdb* files (*polymer_X.pdb* by default). Files are generated in *polymer_molecules* directory

|SMA molecule|DIBMA molecule|
|----|----|
|![SMA](../images/pol_SMA.png)|![DIBMA](../images/pol_DIBMA.png)|

<p align="center">
  <img src="../images/pol_SMA.png" width="300" title="SMA molecule"/>
  <img src="../images/pol_DIBMA.png" width="450" title="DIBMA molecule"/> 
</p>

# Description

*polymer_generation.py* builds polymer molecules consisting of diisobutylene, maleic acid, styrene monomers. The syntax is (default values are displayed):

`polymer_generation.py -n 100 -l 36 -r 3 -p 3 -ps 0 -pm 0.25 -pd 0.75 -o polymer`

Options are:

**-n** Number of polymer molecules to generate

**-meanlen** Mean length of polymer molecules

**-rmsd** Root mean square deviation of length of polymer molecules

**-protstate** Number of protonation state. Protonation state is a number that defines the mean charge of maleic acid monomer. Correspondence of protonation state, pH value and mean charge of MA monomer is given by table:

| Protonation state   | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| pH  | 5 | 6 | 7 | 8 | 9 | 10 | >10 |
| charge | -0.3 | -0.5 | -1.0 | -1.2 | -1.7 | -1.9 | -2.0 |

**-pm** Occurrence of maleic acid monomer

**-pd** Occurrence of diisobutylene monomer

**-ps** Occurrence of styrene monomer

**-out_filename** Name of generated *.pdb* files (*polymer_X.pdb* by default). Files are generated in *polymer_molecules* directory

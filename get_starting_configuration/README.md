# Description

*orient_poly.py* finds principal axes of polymer molecules, alignes the longest axis along z-axis and then places polymer molecules in rings around ridges of lipid membrane. **-n** specifies the number of layers of polymer molecules.
```
python orient_poly.py -m membrane.pdb -p directory_polymer_molecules_reside -o outfile_name.pdb -n 2 --r_coeff 1.0 --angle_coeff 2.0
```
There are two coefficients **--r_coeff** and **--angle_coeff**. *r_coeff = 1.0* and *angle_coeff = 2.0* are the least values when it is garanteed that none of the polymer molecules overlap. Decrease these coefficients for packing polymer molecules closer to each other but always check the conformation obtained because overlapping of the structures may occur.

<p align="center">
  <img width="800" src="../images/starting_configuration.jpg">
</p>

*add_ter_between_chains.py* adds TER strings in *.pdb* files in places where the new count of resid begins.
```
python add_ter_between_chains.py -i infile_name.pdb -o outfile_name.pdb
```

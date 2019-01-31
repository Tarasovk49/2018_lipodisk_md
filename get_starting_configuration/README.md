**Usage**

Script *orient_poly.py* finds principal axes of polymer molecules, alignes the longest axis along z-axis and then places polymer molecules in rings around ridges of lipid membrane. *-n* specifies number of layers.
```
orient_poly.py -m membrane.pdb -p directory_polymer_molecules_reside -o outfile_name.pdb -n 3
```
In the script there are two coefficients *r_coeff* and *angle_coeff*. *r_coeff = 1.0* and *angle_coeff = 2.0* are the least values when it is garanteed that none of the polymers overlap. Decrease those coefficients for packing polymer molecules closer to each other but always check if there is an overlapping of the structures.

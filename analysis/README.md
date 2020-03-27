# Description

This directory provides analysis scripts of Gromacs trajectories for lipodisks with example outputs.

### Calculate area per lipid
Calculation of area per lipid is provided by two different scripts.
* The first one works much faster and calculates area per lipid for projections of lipids onto a fitted to those lipids plane -  [*APL_calculation_planefit.py*](APL_calculation_planefit.py).
* The other one finds projections of lipids onto a fitted to those lipids parabola. That is useful if your membrane bended during the simulation. It is quite normal for membranes with area greater than ~200 nm^2 - [*APL_calculation_parabolafit.py*](APL_calculation_parabolafit.py).
### Calculate order parameters of lipids

### Calculate rotational correlation time of lipids

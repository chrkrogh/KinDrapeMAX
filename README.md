# KinDrapeMAX
Written by Christian Krogh, Ph.D   
Department of Materials and Production, Aalborg University, Denmark.

Licensed under MIT license. Copyright (c) 2023 Christian Krogh

The program can analyze and optimize the draping pattern on a mold using 
a kinematic draping model and the built-in MATLAB implementation of 
genetic algorithm, i.e. the function ga. It is intended for wind turbine 
blade production, i.e. where multiple roll-widths or courses of fabric 
are necessary to cover the mold and where the mold is more or less 
rectangular. It can, however, also analyze a single ply on any smooth 
surface of the form z = F(x,y), i.e. with a unique z coordinate for every coordinate pair (x,y). 
The program is compatible with MATLAB version R2019b through R2022a. It will likely be compatible with future releases as well.
Please see the [user guide](/KinDrapeMAX_user_guide.md) for more information on how to run it.

The program was used to generate the results described in the journal paper Krogh, C., Hermansen, S. M., Lund, E., Kepler, J., & Jakobsen, J. (2023). A matter of course: Generating optimal manufacturing instructions from a structural layup plan of a wind turbine blade. Composites Part A: Applied Science and Manufacturing, 172, 107599. Link: https://doi.org/10.1016/j.compositesa.2023.107599

Below is an example of output generated with the program: five courses of fabric draped in a mold and colored based on shear angle.
![This is an image](/Output/2022-December-19_14-52_1layer_baseline/ShearFig.png)

Feel free to use the program according to the license but remember to give proper credit by citing the project and the paper with DOI numbers
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7525531.svg)](https://doi.org/10.5281/zenodo.7525531)

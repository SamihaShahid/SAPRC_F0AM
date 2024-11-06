# SAPRC_F0AM

This repository contains SAPRC gas phase mechanism (SAPRC-18, SAPRC-22 and [developed furan mechanism](https://pubs.acs.org/doi/abs/10.1021/acsearthspacechem.0c00058)) files converted to F0AM (Framework for 0-D Atmospheric Modeling) input files. The SAPRC box model photolysis calculation FORTRAN codes are rewritten in MATLAB to run with F0AM.


### Citation
* An Update and Evaluation of SAPRC Gas-Phase Mechanism for Phenolic Compounds using Chamber Data and Biomass Burning Plume Simulations, (In preparation). We ask that if you use the python programs to run simulation in F0AM, you cite this paper.

## Organization of the repository

 * [Python_scripts](python_scripts): Contains the Python functions that are used to create the SAPRC mechanism files in F0AM format. 
 * [SAPRC_mechanism](SAPRC_mechanism): Contains SAPRC mechanism in F0AM format and photolysis calculation MATLAB function files.
 * [Simulation_files](Simulation_files): Lagrangian plume simulations.

## Contact
Samiha Binte Shahid (sbint003@ucr.edu)

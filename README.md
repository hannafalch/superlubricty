# superlubricty
There are two folders in this directory:

"lammps_scripts":
- All code used to create simulations for the thesis, along with the dump files from running the simulations can be found in the folder "lammps_scripts". This folder is organised after the parameter that is altered from its default-value, in the simulations. The folder named "Load" thus shows the simulations run with a varying load. For each simulation one can find all dump files, in addition to the files used to create the simulation: granular.in, run.in and myBatch.sh. One general granular.in-file, only containting default values, is also available in the "lammps_scripts"-folder.

"code":
-All general functions used to analyse to code can be found in the folder "code". Here are functions to read the dump-files, calculate the overlap and potential using Hertz contact mechanics, and implementations of the general analytcal solution. 
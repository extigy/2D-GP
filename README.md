# Introduction
2D-GP is a FORTRAN project designed to numerically solve the Gross-Pitaevskii equation (GPE) in two dimensions (2D). Solving the GPE allows for qualitatively accurate simulations of Bose-Einstein Condensates (BECs) at zero temperature.  
2D-GP solves the GPE using 4th order Runge-Kutta time stepping on a grid of points with regular spacing in both the x and y dimensions. User defined grid spacing and time step is supported.

### Requirements
2D-GP requires NetCDF and NetCDF-Fortran installations to run. NetCDF is used to save compressed data files. On Ubuntu you can install the required packages by running:
```
sudo apt-get install libnetcdf11 libnetcdff6 libnetcdf-dev libnetcdff-dev
```

You can also create a local installation of NetCDF by downloading the source files from the [NetCDF website](https://www.unidata.ucar.edu/software/netcdf/) and compiling them yourself. If you do this, make a note of the installation location as you will need it to install 2D-GP.

### Installation
* Clone the git repository somewhere : `git clone https://github.com/Extigy/2D-GP.git`

* Run `./install` to setup. You can also run  `./install <install-dir>` to install to any desired installation location.

* If you have a non-system installation of NetCDF (i.e. you compiled it from source), edit the file `Makefile` and at the top change the `NETCDF` and `NETCDF-FORTRAN` variables to match your installation directory.

* New terminals should now be able to run  `make2dgp` anywhere.

### Running a Simulation
* Create a new simulation directory
* Enter the directory and run `make2dgp` to set up a simulation at your current location.
* Type `./gp&` to start the simulation.
* The simulation status is printed to the file `STATUS`.

### Editing Parameters
To run a simulation with custom parameters
*  Create a simulation as in **Running a Simulation** but do not run `./gp&`.
*  Edit `params.in` to include any parameters you wish to change from their defaults and run `make2dgp` again to reflect the changes.
*  Run `./gp&` to start the simulation.

A full list of parameters, their descriptions, and their default values is included at the end. The code is written so that, most editable parameters can be left as default and so not included in `params.in` at all.

# Documentation
See [DOCUMENTATION.html](https://rawgit.com/extigy/2D-GP/master/DOCUMENTATION.html).
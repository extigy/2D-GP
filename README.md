# Introduction
2D-GP is a FORTRAN project designed to numerically solve the Gross-Pitaevskii equation (GPE) in two dimensions (2D). Solving the GPE allows for qualitatively accurate simulations of Bose-Eintein Condensates (BECs) at zero temperature.

### Numerical Methods
2D-GP solves the GPE using 4th order Runge-Kutta time stepping on a grid of points with regular spacing in both the x and y dimensions. User defined grid spacing and time step is supported.

### Installation
* Clone the git repository: `git clone https://github.com/Extigy/2D-GP.git`

* Run `./install` to setup. You can also run  `./install <install-dir>` to install to any desired installation location.

* New terminals should now be able to run  `./make2dgp` anywhere.

### Running a Simulation
* Create a new simulation directory
* Enter the directory and run `make2dgp` to set up a simulation at your current location.
* Type `./gp&` to start the simulation.

### Editing Parameters
To run a simulation with custom parameters
*  Edit `params.in` and run `make2dgp` again to reflect the changes.
*  Run `./gp&` to start the simulation.

A full list of editable parameters can be found at the end of this document. Hopefully, most editable parameters can be left as default.

# Gross-Pitaevskii Equation
2D-GP solves the GPE in 2 distinct and different dimensionless forms: the homogenous system and the harmonic trapped system.  

The homogeneous system is used to fill the entire computational box with fluid, emulating superfluid liquid helium II. The harmonic trapped case is employed to emulate quasi-2d BEC experiments when confined using a harmonic trap.

Both forms of the governing equation are shown below.
### Homogenoeus System
The dimensionless homogenous GPE is defined as   
![H_GPE](http://latex.codecogs.com/gif.latex?i\frac{\partial%20\psi}{\partial%20t}%20%3D%20-\frac{1}{2}%20\nabla^2\psi%20&plus;%20|\psi|^2\psi%20&plus;%20V\psi%20-\psi + iv_{ob}\frac{\partial%20\psi}{\partial%20x}),  
where ![psi](http://raw.githubusercontent.com/Extigy/2D-GP/params/images/psi.gif) is the condensate wavefunction, *V* is a potential allowing for obstacles and ![v_ob](http://github.com/Extigy/2D-GP/blob/params/images/v_ob.gif) is the fluid velocity along the x direction.  

The above version of the GPE is valid when properties are scaled with the so called `natural units',

![Natural Scaling](http://raw.githubusercontent.com/Extigy/2D-GP/params/images/natural_scale.gif)
-----
### Harmonically Trapped System
The dimensionless GPE in a harmonically trapped system is defined as  
![T_GPE](http://latex.codecogs.com/gif.latex?i\frac{\partial%20\psi}{\partial%20t}%20%3D%20-\frac{1}{2}%20\nabla^2\psi%20&plus;%20g_{2D}|\psi|^2\psi%20&plus;%20V\psi%20-\mu_{2D}\psi),  
where ![psi](http://raw.githubusercontent.com/Extigy/2D-GP/params/images/psi.gif) is the condensate wavefunction,...

The above version of the GPE is valid when...

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
2D-GP solves the GPE in 2 different dimensionless forms: the homogeneous system and the harmonically trapped system.  

The homogeneous system is used to fill the entire computational box with fluid, emulating superfluid liquid helium II. The harmonic trapped case is employed to emulate quasi-2d BEC experiments when confined using a harmonic trap.

---
### Homogeneous System
The dimensionless homogeneous GPE is defined as
![H_GPE](http://latex.codecogs.com/gif.latex?i\frac{\partial\psi}{\partial%20t}-\frac{1}{2}\nabla^2\psi+|\psi|^2\psi+V\psi-\psi+iv_{ob}\frac{\partial\psi}{\partial%20x}),  
where ![psi](http://latex.codecogs.com/gif.latex?%5Cpsi) is the condensate wavefunction, ![V](http://latex.codecogs.com/gif.latex?V) is a potential allowing for obstacles and ![v_ob](http://latex.codecogs.com/gif.latex?v_{ob}) is the fluid velocity along the x direction.  

The above version of the GPE is valid when properties are scaled with the so called '**natural units**',  
![Natural Scaling](http://latex.codecogs.com/gif.latex?%5Cdpi%7B110%7D%20%5C%5C%5Cmathrm%7BDensity%7Eat%7Einfinity%3A%7D%7En_0%20%5C%5C%5Cmathrm%7BLength%3A%7D%7E%5Cxi%20%3D%20%5Cfrac%7B%5Chbar%7D%7B%5Csqrt%7Bm%5Cmu%7D%7D%20%5C%5C%5Cmathrm%7BEnergy%3A%7D%7E%5Cmu%20%3D%20n_0g%20%5C%5C%5Cmathrm%7BVelocity%3A%7D%7Ec%3D%5Cfrac%7B%5Csqrt%7Bn_0g%7D%7D%7Bm%7D%20%5C%5C%5Cmathrm%7BTime%3A%7D%7E%5Cfrac%7B%5Cxi%7D%7Bc%7D)

The following parameters are related to the homogeneous equation.

Parameter | Default Value | Explanation
--- | --- | ---
`RHSType` | `0` | GPE Type - Set to `0` for **natural units**, `1` for **harmonic oscillator units**.
`VOBS` | `0` | First simulation's velocity. ![v_ob](http://latex.codecogs.com/gif.latex?v_{ob}=) ![v_ob](http://latex.codecogs.com/gif.latex?\dpi{80}\mathrm{VOBS}/100).
`VOBE` | `0` | Final simulation's velocity. This is ![v_ob](http://latex.codecogs.com/gif.latex?v_{ob}=) ![v_ob](http://latex.codecogs.com/gif.latex?\dpi{80}\mathrm{VOBE}/100).
`VOBST` | `1` | Increase ![v_ob](http://latex.codecogs.com/gif.latex?v_{ob}) by ![v_ob](http://latex.codecogs.com/gif.latex?\dpi{80}\mathrm{VOBST}/100) per simulation.
***
### Harmonically Trapped System
The dimensionless GPE in a harmonically trapped system is defined as  
![T_GPE](http://latex.codecogs.com/gif.latex?i\frac{\partial%20\psi}{\partial%20t}%20%3D%20-\frac{1}{2}%20\nabla^2\psi%20&plus;%20g_{2D}|\psi|^2\psi%20&plus;%20V\psi%20-\mu_{2D}\psi),  
where ![psi](http://latex.codecogs.com/gif.latex?\psi) is the condensate wavefunction,...

The above version of the GPE is valid when...

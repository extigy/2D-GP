# Introduction
2D-GP is a FORTRAN project designed to numerically solve the Gross-Pitaevskii equation (GPE) in two dimensions (2D). Solving the GPE allows for qualitatively accurate simulations of Bose-Eintein Condensates (BECs) at zero temperature.  
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
*  Create a simulation as in **Running a Simulation** but do not run `./gp&`.
*  Edit `params.in` and run `make2dgp` again to reflect the changes.
*  Run `./gp&` to start the simulation.

Editable parameters along with descriptions and their default values can be found throughout this document. Hopefully, most editable parameters can be left as default.

# Gross-Pitaevskii Equation
2D-GP solves the GPE in 2 different dimensionless forms: the homogeneous system and the harmonically trapped system.  

The homogeneous system is used to fill the entire computational box with fluid, emulating superfluid liquid helium II. The harmonic trapped case is employed to emulate quasi-2D BEC experiments when confined using a harmonic trap.

---
### Homogeneous System
The dimensionless homogeneous GPE is defined as
![H_GPE](http://latex.codecogs.com/gif.latex?i%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%20t%7D=-%5Cfrac%7B1%7D%7B2%7D%5Cnabla%5E2%5Cpsi%2B%7C%5Cpsi%7C%5E2%5Cpsi%2BV%5Cpsi-%5Cpsi%2Biv_%7Bob%7D%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%20x%7D),  
where ![psi](http://latex.codecogs.com/gif.latex?%5Cpsi) is the condensate wavefunction, ![V](http://latex.codecogs.com/gif.latex?V) is a potential allowing for obstacles and ![v_ob](http://latex.codecogs.com/gif.latex?v_{ob}) is the fluid velocity along the x direction.  

This version of the GPE is valid when properties are scaled with the so called '**natural units**'


![Natural Scaling](http://latex.codecogs.com/gif.latex?%5Cdpi%7B110%7D%20%5C%5C%5Cmathrm%7BDensity%7Eat%7Einfinity%3A%7D%7En_0%20%5C%5C%5Cmathrm%7BLength%3A%7D%7E%5Cxi%20%3D%20%5Cfrac%7B%5Chbar%7D%7B%5Csqrt%7Bm%5Cmu%7D%7D%20%5C%5C%5Cmathrm%7BEnergy%3A%7D%7E%5Cmu%20%3D%20n_0g%20%5C%5C%5Cmathrm%7BVelocity%3A%7D%7Ec%3D%5Cfrac%7B%5Csqrt%7Bn_0g%7D%7D%7Bm%7D%20%5C%5C%5Cmathrm%7BTime%3A%7D%7E%5Cfrac%7B%5Cxi%7D%7Bc%7D)

Running multiple simulations sequentially at varying velocities is supported. The solver will first run a simulation at the ![v_ob](http://latex.codecogs.com/gif.latex?v_{ob}) defined by `VOBS`.
Once complete the solver will increase ![v_ob](http://latex.codecogs.com/gif.latex?v_{ob}) by the amount defined by `VOBST` and rerun the simulation. This will continue until ![v_ob](http://latex.codecogs.com/gif.latex?v_{ob})
is larger than the value defined by `VOBE` at which point the solver will terminate.

The following parameters can be modified,

Parameter | Default | Explanation
--- | --- | ---
`RHSType` | `0` | GPE Type - Set to `0` for **natural units**, `1` for **harmonic oscillator units**.
`VOBS` | `0` | First simulation's velocity. ![v_ob](http://latex.codecogs.com/gif.latex?v_{ob}=) ![v_ob](http://latex.codecogs.com/gif.latex?%5Cdpi{80}~%5Cmathrm{VOBS}/100).
`VOBE` | `0` | Final simulation's velocity. This is ![v_ob](http://latex.codecogs.com/gif.latex?v_{ob}=) ![v_ob](http://latex.codecogs.com/gif.latex?%5Cdpi{80}~%5Cmathrm{VOBE}/100).
`VOBST` | `1` | Increase ![v_ob](http://latex.codecogs.com/gif.latex?v_{ob}) by ![v_ob](http://latex.codecogs.com/gif.latex?%5Cdpi{80}%5Cmathrm{VOBST}/100) per simulation.
***
### Harmonically Trapped System
The dimensionless GPE in a harmonically trapped system is defined as  
![T_GPE](http://latex.codecogs.com/gif.latex?i%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%20t%7D=-%5Cfrac%7B1%7D%7B2%7D%5Cnabla%5E2%5Cpsi%2Bg_{2D}%7C%5Cpsi%7C%5E2%5Cpsi%2BV%5Cpsi-%5Cmu_{2D}%5Cpsi),  
where ![psi](http://latex.codecogs.com/gif.latex?%5Cpsi) is the condensate wavefunction, ![g2d](http://latex.codecogs.com/gif.latex?g_{2D}) is the interaction strength, ![V](http://latex.codecogs.com/gif.latex?V) is a potential allowing for a harmonic trap and/or obstacles, and ![mu2d](http://latex.codecogs.com/gif.latex?%5Cmu_{2D}) is the chemical potential.

This version of the GPE is valid when properties are scaled with the so called '**harmonic oscillator units**'


![HO_Scaling](http://latex.codecogs.com/gif.latex?%5Cdpi{120}%20%5C%5C%5Cmathrm{Energy%3A}~%5Cmu%20%3D%20%5Chbar%20%5Comega_r%20%5C%5C%5Cmathrm{Length%3A}~%5Cl_r%20%3D%20%5Csqrt{%5Cfrac{%5Chbar}{m%5Comega_r}}%20%5C%5C%5Cmathrm{Time%3A}~%20%5Comega_r^{-1},)

where ![omega_r](http://latex.codecogs.com/gif.latex?%5Comega_r) is the trap frequency in the radial direction.  
*NB. Strictly speaking, a 2D mean-field description of an oblate trapped condensate is only valid when ![condition](http://latex.codecogs.com/gif.latex?%5Cdpi{100}Nal_z^3/l_r^3%5Cll1), where N is the number of atoms, a is the s-wave scattering length and ![l_z](http://latex.codecogs.com/gif.latex?%5Cdpi{100}%5Cl_z%20%3D%20%5Csqrt{%5Chbar/m%5Comega_z}) and ![l_r](http://latex.codecogs.com/gif.latex?%5Cdpi{100}%5Cl_r%20%3D%20%5Csqrt{%5Chbar/m%5Comega_r}) are the axial and radial harmonic oscillator lengths.*

The following parameters are related to the harmonically trapped system.

Parameter | Default | Explanation
--- | --- | ---
`RHSType` | `0` | GPE Type - Set to `0` for **natural units**, `1` for **harmonic oscillator units**.
`harm_osc_C` | `100` | Value of ![g2d](http://latex.codecogs.com/gif.latex?g_{2D})
`harm_osc_mu` | `5` | Value of ![mu2d](http://latex.codecogs.com/gif.latex?%5Cmu_{2D})
`enableTrap` |`.true.`| Enable or disable the trapping potential
---
### Damped GPE
For both cases of the GPE, an optional damping parameter can be applied, allowing for simulations using the Damped GPE. The damped version of the GPE as a rough and rudimentary simulation of finite temperature effects on BECs.  
The damped GPE is found by replacing the left hand side of the GPE with,  
![T_GPE][dgpe]
[dgpe]: http://latex.codecogs.com/gif.latex?(i-%5Cgamma)%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%20t%7D=
where ![gamma](http://latex.codecogs.com/gif.latex?%5Cgamma) is a small positive real parameter controlling the strength of the damping.

In damped GPE simulations, the condensate norm can also be renormalised at every iteration, preventing atom loss for large damping parameters.  

Parameter | Default | Explanation
--- | --- | ---
`GAMMAC` | `0` | Value of ![gamma](http://latex.codecogs.com/gif.latex?%5Cgamma)
`rtNorm` |`.false.`| If `.true.` the wavefunction is normalised at every time step
---
# Numerical Method
2D-GP uses 4th order Runge-Kutta time stepping on a grid of points with regular spacing in both the x and y dimensions.  
The solver begins by finding an inital condition using the *imaginary time method*, which obtains a steady solution known as the *ground state*. An optional random perturbation can then be applied to this initial condition.  
The solver then runs in *real time* producing a solution to the GPE for the given parameters.


The following parameters can be customised,  

Parameter | Default | Explanation
--- | --- | ---
`NX` | `512` | Number of grid points in the x direction
`NY` | `512` | Number of grid points in the y direction
`DSPACE` | `0.05` | Grid spacing in dimensionless units
`DTSIZE` | `0.001` | Time step size in dimensionless units
`ISTEPS` | `2000` | Number of iterations to run the solver in the *imaginary time method*
`NSTEPS` | `10000` | Number of iterations to run the solver in real time
`noiseamp` | `0` | Amplitude of random noise to add to the ground state initial condition

### Boundary Conditions
Both **reflective** and **periodic** boundary conditions are supported.

Parameter | Default | Explanation
--- | --- | ---
`BCX` | `0` | Boundary conditions in the x direction - Set to `0` for **reflective**, `1` for **periodic**.
`BCY` | `0` | Boundary conditions in the y direction - Set to `0` for **reflective**, `1` for **periodic**.

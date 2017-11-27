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

# Gross-Pitaevskii Equation
2D-GP solves the GPE in 2 different dimensionless forms: the homogeneous system and the harmonically trapped system.  

The homogeneous system is primarily used to emulate superfluid helium II. The harmonic trapped case is primarily used to emulate quasi-2D BEC experiments, confined using a harmonic trap.

### Homogeneous System
The dimensionless homogeneous GPE is defined as
$$i\frac{\partial\psi}{\partial t}=-\frac{1}{2}\nabla^2\psi+|\psi|^2\psi+V\psi-\psi+iv_{ob}\frac{\partial\psi}{\partial x},$$
where $\psi$ is the wavefunction, $V$ is a potential allowing for the use of boundaries and obstacles and $v_{ob}$ is the fluid velocity along the $x$ direction.  

In this setup the wavefunction is scaled so that the density away from any perturbations is $n_0=1$.

This version of the GPE is valid when properties are scaled with the so called '**natural units**',
 * Density at infinity: $n_0$,
 * Length: $\xi = \dfrac{\hbar}{\sqrt{m\mu}}$,
 * Energy: $\mu = n_0g$,
 * Velocity: $c = \dfrac{\sqrt{n_0g}}{m}$,
 * Time: $\dfrac{\xi}{c}$.

To use this system, set `RHSType=0` in `params.in`.

### Harmonically Trapped System
The dimensionless GPE in a harmonically trapped system is defined as 
$$i\frac{\partial\psi}{\partial t}=-\frac{1}{2}\nabla^2\psi+g_{2D}|\psi|^2\psi+V\psi-\mu_{2D}\psi - i \hbar\Omega\(x\frac{\partial \psi}{\partial y}-y\frac{\partial \psi}{\partial x}\),$$
where $\psi$ is the condensate wavefunction, $g_{2D}$ is the interaction strength, $V$ is a potential allowing for the use of boundaries and obstacles, $\mu_{2D}$ is the chemical potential, and $\Omega$ is the angular momentum of the rotating frame.

This version of the GPE is valid when properties are scaled with the so called '**harmonic oscillator units**':
 * Length: $\displaystyle l_r = \sqrt{\frac{\hbar}{m\omega_r}}$,
 * Energy: $\mu = \hbar\omega_r$,
 * Time: $\omega_r^{-1}$,

where $\omega_r$ is the trap frequency in the radial direction. To use this system, set `RHSType=1` in `params.in`.

***Note**: Strictly speaking, a 2D mean-field description of an oblate trapped condensate is only valid when $\dfrac{Nal_z^3}{l_r^3}\ll1$, where N is the number of atoms, $a$ is the s-wave scattering length and $l_z=\sqrt{\dfrac{\hbar}{m\omega_z}}$ and $l_r = \sqrt{\dfrac{\hbar}{m\omega_r}}$ are the axial and radial harmonic oscillator lengths.*


### Damped GPE
In both the harmonic trap system and the homogeneous system, an optional damping parameter can be applied. The damped version of the GPE is a very rough phenomenological simulation of finite temperature effects on BECs.

The damped GPE replaces the left hand side of the GPE with $(i-\gamma)\dfrac{\partial\psi}{\partial t}$, where $\gamma$ is a small positive parameter controlling the rate of damping.

The value of $\gamma$ can be set to 0.01, for example, by adding `GAMMAC=0.01d0` to `params.in`.

# The Numerical Method
To solve the GPE, 2D-GP uses 4th order Runge-Kutta time stepping on a grid of points with regular spacing in both the $x$ and $y$ dimensions. There are three main modes of time stepping. `ISTEPS` controls how many steps to run *imaginary time propagation*, `DSTEPS` controls how many steps to run with a *highly damped propagation*, and `NSTEPS` controls how many steps of *real time propagation* to run.

You will want `ISTEPS` and `DSTEPS` to be small or zero, used simply to setup the system, with the majority of time spent running `NSTEPS` of real-time simulation.

### Imaginary Time Propagation

For imaginary time propagation we move from real to imaginary time using the substitution $t_i = -it$. This transforms the GPE into a form similar to a diffusion equation. As a result, a local equilibrium can be found by propagating through time. This causes the system to approach a ground state solution. The cost of using a diffusion equation is that solution decays over time. We inhibit the overall decay of the wavefunction by renormalising during imaginary time propagation. The imaginary time propagation method converges on the ground state solution very slowly and so we must perform many numerical steps to prepare an initial state.

The renormalisation step forces the solution to have a constant norm throughout imaginary time propagation. The norm that the system is renormalised to depends on the choice of GPE units, but it is fundamentally an approximation based on the size of the condensate and the imposed potential, $V$. This works well when $\gamma=0$, because then the energy and norm are conserved during real time simulation. However, when $\gamma \ne 0$ the system is slightly dissipative and so the energy and norm change over time to approach the energy and norm of the true ground state solution. If the norm approximation by 2D-GP is particularly bad, this can be noticeable.

Imaginary time propagation is run before the initial condition defined in the file `ic.in` is multiplied into the system.

### Highly Damped Propagation

After imaginary time propagation, initial conditions (defined in `ic.in`) are applied and the system is run with a highly damped GPE for `DSTEPS` steps. During this period the system will continue to approach the ground state. Running highly damped propagation for a short amount of time has two benefits:
 * If your initial condition in `ic.in` is a combination of several solutions multiplied into the system, the resulting wavefunction is most likely not a true solution of the GPE and so an immediate emission of sound can occur. With highly damped propagation the initial start to the simulation is 'smoothed' as any emitted sound is damped away
 * If the norm approximation used by 2D-GP in imaginary time propagation is particularly bad, running the system highly damped for a short amount of time forced the norm of the system to approach the norm of the true ground state solution.

By default, the phase of the system is continually enforced during this period. This can be disabled by setting `HOLDIC = .false.`. However, in this case `DSTEPS` should be small if you are using a complicated `ic.in`. Vortices, in particular, will approach each other and annihilate unless the phase of the system is enforced throughout.

### Real Time Propagation

After highly damped propagation, numerical noise is optionally applied to break symmetry (with e.g. `noiseamp=0.05d0` in `params.in`) and the system is run normally for `NSTEPS` steps. This is where the majority of GPE simulation occurs and by default it is only the output wavefunction data from this mode that is written to disk.

### Boundary Conditions
**Reflective**, **periodic** and **zero** boundary conditions are supported and can be independently applied to the $x$ and $y$ directions.

`BCX` controls the boundary conditions in the $x$ direction and `BCY` controls the boundary conditions in the $y$ direction. Set to 0 for reflective, 1 for periodic, 2 for zero.

# Adding an Initial Condition
An initial condition can be imprinted (such as imprinting vortices) at the start of the simulation. `ISTEPS` and `DSTEPS` should be tweaked to reflect how much imaginary time propagation you would like before and after the imprinting.

To add an initial condition, after running `make2dgp` in an empty directory edit the file `ic.in` with Fortran code to be run. Some vortex imprinting functions are provided in the file `utils.f95`. An example ic.in could be:
```
!Enter initial conditions in the form
!call IC(...)
!
call insert_vortex(-10.0d0,0.0d0,1)
call insert_vortex(10.0d0,0.0d0,-1)
```
This imprints a vortex phase and density profile at $x,y=(-10,0)$ with positive circulation, and a vortex phase and density profile at $x,y=(10,0)$ with negative circulation.

The `GRID` Fortran variable contains the system wavefunction, it is this variable you should make changes to apply some initial condition. For a complex setup involving temporary variables, you should create a new function in `utils.f95` (using `insert_vortex` as a template) and call your new function from `ic.in`.

# Trapping and Obstacle Potentials
Trapping and obstacle potentials are controlled separately and combined by 2D-GP to form the potential term in the GPE, $V = V_{\text{trap}}+V_{\text{obj}}$. The potential is measured in units of $\mu$.

## Trapping Potential $V_{\text{trap}}$

To enable a trapping potential, set `enableTrap=.true.` in `params.in`. The type of trapping potential is chosen by setting the `trapType` variable.

### Harmonic Trap
`trapType=0` sets up a standard trapping potential for use with a harmonically trapped system, $$V_{\text{trap}} = \frac{r^2}{2}.$$

### Hard Circle Trap
`trapType=1` sets up a hard wall potential that takes the value of 0 inside the trap boundary and some maximum value outside. The shape of the trap boundary is a circle,
$$V_{\text{trap}} = \begin{cases}
T_{\text{max}},  & \text{if $|r|>T_r$,}\cr
0, & \text{otherwise.}
\end{cases}$$
where $T_{\text{max}}$ (`TRAPHEIGHT`) is the maximum value of the trapping potential and $T_r$ (`TRAPR`) is the radius of the boundary.

### Hard Box Trap
`trapType=2` sets up a hard wall potential that takes the value of 0 inside the trap boundary and some maximum value outside. The shape of the trap boundary is a box,
$$V_{\text{trap}} = \begin{cases}
T_{\text{max}},  & \text{if $|x|>T_r$ or $|y|>T_r$,}\cr
0, & \text{otherwise,}
\end{cases}$$
where $T_{\text{max}}$ (`TRAPHEIGHT`) is the maximum value of the trapping potential and $T_r$ (`TRAPR`) is a the size of the bounding box.

### Soft Circle Trap
`trapType=3` sets up a smoothly varying potential that takes the value of 0 inside the trap boundary and some maximum value outside. The shape of the trap boundary is a circle,
$$V_{\text{trap}} = \frac{T_{\text{max}}}{1+e^{T_{\beta}(T_r-|r|)}},$$
where $T_{\text{max}}$ (`TRAPHEIGHT`) is the maximum value of the trapping potential, $T_r$ (`TRAPR`) is the radius of the bounding circle, and $T_{\beta}$ (`TRAPBETA`) controls the steepness of the soft wall.

### Soft Box Trap
`trapType=4` sets up a smoothly varying potential that takes the value of 0 inside the trap boundary and some maximum value outside. The potential is defined in terms of 
$$V_{\text{box}} = \frac{T_{\text{max}}}{1+e^{T_{\beta}(T_r-|x|)}} + \frac{T_{\text{max}}}{1+e^{T_{\beta}(T_r-|y|)}}.$$
So that the shape of the trap boundary is a box,
$$V_{\text{trap}} = \begin{cases}
T_{\text{max}},  & \text{if $V_{\text{box}}>T_{\text{max}}$}\cr
V_{\text{box}}, & \text{otherwise,}
\end{cases}$$
where $T_{\text{max}}$ (`TRAPHEIGHT`) is the maximum value of the trapping potential, $T_r$ (`TRAPR`) is the radius of the bounding circle, and $T_{\beta}$ (`TRAPBETA`) controls the steepness of the soft wall.


### The Shin Experiment

An optional mode can be enabled emulating an experiment where an oblate BEC was dragged around a stationary laser obstacle by translating the trap.
The simulation starts with a stationary trap. The velocity of the trap is increased to its maximum, and finally is slowed down back to stationary after a fixed period of time.

Parameter  | Explanation
--- | ---
`doShin` | Set to `.true.` to enable the Shin experiment.
`TVXDASH` | Trap's maximum velocity in the *x* direction.
`TVYDASH` | Trap's maximum velocity in the *y* direction.
`TTM` | Amount of time the trap will be moving.

## Object Potential $V_{\text{obj}}$

To enable an object potential, set `enablePot=.true.` in `params.in`. The type of trapping potential is chosen by setting the `potType` variable.

### Gaussian Object
`potType=0` enables a Gaussian Obstacle of the form
$$V_{\text{obj}} = V_{\text{max}} \exp \left[ \frac{(x-x_0)^2}{r_x^2} + \frac{(y-y_0)^2}{r_y^2} \right],$$
where $r_x$ (`RRX`) and $r_y$ (`RRY`) are the obstacle radius along the $x$ and $y$ directions and $(x_0,y_0)$ (`OBJXDASH`, `OBJYDASH`) is the obstacle position. $V_{\text{max}}$ (`OBJHEIGHT`) is the maximum height of the potential.

#### Oscillation
The Gaussian obstacle can forced to oscillate, emulating wires or tuning forks as used in liquid helium experiments. The oscillation is achieved by a coordinate transform for the obstacle only of $x\rightarrow x+A\sin(\omega t)$, where $A$ (`OBJAMP`) is the oscillation amplitude and $\omega$ (`OBJW`) the frequency.

#### Rotating
`potType=1` enables a rotating Gaussian obstacle. The obstacle can rotate freely about its centre, and reacts to forces and pressures exerted by the fluid.  Here the obstacle has an initial angle (`OBJANGLE`), initial angular velocity (`OBJANGLEV`) and a [moment of inertia](http://en.wikipedia.org/wiki/Moment_of_inertia) (`MOMINERTIA`).

### Bitmap Sourced Potential
`potType=4` enables the loading of a bitmap image which is then used to define a potential. To use this feature first create a grayscale bitmap of size (`NX+1` $\times$ `NY+1`). White areas correspond to areas where the potential will be equal to $V_{\text{max}}$ (`OBJHEIGHT`), while black areas correspond to areas where the potential will be 0. Grayscale colours will cause the potential to scale between 0 and $V_{\text{max}}$.

Set the parameter `pot_filename` equal to the location of the potential `.bmp` image.

# Output Data
Data is output every set amount of time steps. This output frequency can be easily customised. `dumpwf` sets the wavefunction output frequency and `dumputil` sets the utility output frequency.

### Wavefunction Output

Wavefunction files are output with the filename: `dumpwf.ssssss.nc`. Here `ssssss` is the number of the output file. e.g `dumpwf.000000.nc`.

These files are formatted as NetCDF files and contains the following data:

Name | Size |Description
--- | --- |--- 
`x` | `NX+1` |Grid points along the x dimension
`y` | `NY+1`  |Grid points along the y dimension
`real` | (`NX+1` $\times$ `NY+1`) |Real part of the system wavefunction
`imag` | (`NX+1` $\times$ `NY+1`) |Imaginary part of the system wavefunction
`pot` | (`NX+1` $\times$ `NY+1`) |The current external potential field.

#### Reading Wavefunction Output
Included in the `scripts` folder is a matlab script named `gpeget2dPSI.m`. Use this file to load the output wavefunction data into matlab for post processing. The function takes the simulation directory and the number of the data file.
```
[gridx,gridy,psi,potential] = gpeget2dPSI('/directory/to/simulation',250);
```

### Utility Output

The energy and norm are calculated and output by default every 100 time steps as diagnostic check as the simulation is running. 2D-GP will create a file at the start of the simulation named `utils.dat`, and continuously write to the file with:

Column 1 | Column 2 | Column 3 | Column 4
--- | --- | --- | ---
$-it$, imaginary time | $t$, simulation time | $E$, the total energy at time $t$ | $N$, the norm at time $t$

# Example `params.in`

An example (and default) `params.in` is provided by 2D-GP. It sets up a trapped BEC system and and looks like this:
```
!Enter custom parameters in the form
!PARAMETER_NAME = VALUE
!
NX = 256
NY = 256
NSTEPS=5000
DSTEPS=1000
ISTEPS=100

GAMMAC = 0.1d0
DSPACE = 0.1d0
DTSIZE = 0.002d0

RHSType = 1
BCX = 1
BCY = 1
harm_osc_C = 2000.0d0
harm_osc_mu = 25.26698674d0

enableTrap = .true.
trapType = 0
```

## List of Parameters

Parameter | Default, if omitted | Explanation
--- | --- | ---
`NSTEPS` | `1000` | Number of iterations to run the solver in real time.
`ISTEPS` | `0` | Number of iterations to run the solver in imaginary time.
`DSTEPS` | `0` | Number of iterations to run the solver highly damped.
`HOLDIC` | `.true.` | If true, the phase defined in `ic.in` is constantly re-imprinted during `ISTEPS` and `DSTEPS`.
`PLOTALL` | `.false.` | If true, wavefunction files are written during `ISTEPS` and `DSTEPS`.
`NX` | `512` | Number of grid points in the x direction, not including an extra grid point for zero.
`NY` | `512` | Number of grid points in the y direction, not including an extra grid point for zero.
`DSPACE` | `0.1` | Grid spacing in dimensionless units
`DTSIZE` | `0.002` | Time step size in dimensionless units
`dumpwf` | `100` | Wavefunction output frequency
`dumputil` | `100` | Force/energy output frequency
`RHSType` | `1` | GPE Type - `0` for **natural units**, `1` for **harmonic oscillator units**.
`harm_osc_C` | `2000` | Value of $g_{2D}$.
`harm_osc_mu` | `25.267` | Value of $\mu_{2D}$.
`GAMMAC` | `0` | Value of $\gamma$.
`NORMALL` | `.false` | If true, the GPE is renormalised every time step.
`BCX` | `0` | Boundary conditions in the x direction.
`BCY` | `0` | Boundary conditions in the y direction.
`noiseamp` | `0` | Amplitude of random noise to add to the initial condition.
`ROM` | `0` | Angular momentum of the rotating frame, $\Omega$
`VOB` | `0` | Velocity of the linearly translating frame, $v_{ob}$
`enablePot` | `.true.` | Enable the potential obstacle.
`potType` | `-1` | Obstacle type.
`RRX`|`2.0`| Value for $r_x$, obstacle radius in the $x$ direction.
`RRY`|`2.0`| Value for $r_y$, obstacle radius in the $y$ direction.
`OBJHEIGHT`| `50` | Value for $V_{\text{max}}$, the maximum barrier height. A maximum height of 0 is equivalent to having no obstacle.
`OBJXDASH` | `0` | Value for $x_0$, the obstacle's $x$ position.
`OBJYDASH` | `0` | Value for $y_0$, the obstacle's $y$ position.
`OBJAMP`| `0` | Oscillating obstacle's amplitude.
`OBJW`| `0` | Oscillating obstacle's frequency.
`OBJANGLE`|`0.1`|Initial obstacle angle for rotating obstacle.
`OBJANGLEV`|`0.1`|Initial angular velocity for rotating obstacle.
`MOMINERTIA`|`0.000001`| Rotating obstacle's [Moment of Inertia](http://en.wikipedia.org/wiki/Moment_of_inertia)
`enableTrap` | `.true.` | Enable the potential trap.
`trapType` | `-1` | Trap type.
`TRAPHEIGHT`| `50` | Value for $T_{\text{max}}$, the maximum trap height. A maximum height of 0 is equivalent to having no trapping potential.
`TXDASH`|`0`| Value for $x_0$, center of the trap in $x$.
`TYDASH`|`0`| Value for $y_0$, center of the trap in $y$.
`TXSCALE`| `1` | A scaling factor for the trap in $x$.
`TYSCALE`| `1` | A scaling factor for the trap in $y$.
`TRAPR`| `5` | Box/Circle trap size.
`TRAPBETA`| `2` | Soft wall trap steepness.
`potRep` | `.false.` | If true, recalculates the potential term $V$ at every time step.
`doShin` | `.false.` | If true, runs a simulation of the Shin experiment.
`TTM`| `0` | Amount of time the trap will be moving in the Shin experiment.
`TVXDASH` | `0` | Trap's maximum velocity in the $x$ direction in the Shin experiment.
`TVYDASH`| `0` | Trap's maximum velocity in the $y$ direction in the Shin experiment.
`pot_filename` | '' | Filename for bitmap sourced potential.


# TODO: Vortex Killing and damped region
```
  !Vortex Killer
  logical :: doVortexKilling = .false.
  logical :: vortexKillX = .true.
  logical :: vortexKillY = .false.
  double precision :: vortexKillXDist = 40.0d0
  double precision :: vortexKillYDist = 0.0d0
  integer :: vortexKillFreq = huge(100)
  integer :: killgamma = 0
  !Damping radius
  logical :: dampedX = .false.
  logical :: dampedY = .false.
  logical :: dampedR = .false.
  double precision :: dampedXDist = 10.0d0
  double precision :: dampedYDist = 0.0d0
  double precision :: dampedRDist = 0.0d0
  double precision :: dampedRWidth = 10.0d0
  double precision :: dampedGamma = 0.0d0
```

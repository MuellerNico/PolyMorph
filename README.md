# PolyMorph

TODO

## Abstract

PolyMorph is an extension of PolyHoop, a 2D mechanical tissue simulation software for soft particle dynamics developed by Vetter et al. PolyMorph adds a numerical solver for reaction-diffusion equations on top of it which allows, among other things, the modeling of morphogenetic problems with realistic cell shapes. 

The simulation software is programmed in C++. It makes use of finite-difference approximations, explicit Euler time integration, scatter-gather interpolation between Eulerian and Lagrangian reference frames and shared memory parallelism with OpenMP. 

The solver's correctness and convergence was verified and its computational performance quantified through parallel scaling analysis and runtime comparison to the purely mechanical simulation. The coupled system is then used to study patterning precision with noisy morphogen gradients and the results are compared to previous findings.

## Usage Instructions: Quick Start 
```shell
$ git clone https://github.com/MuellerNico/PolyMorph.git
$ cd PolyMorph/
$ make
$ OMP_NUM_THREADS=8 ./polymorph
```

## Features
- Most PolyHoop capabilities (see limitations)
- Finite difference solver for reaction-diffusion equations
- Advection and dilution terms (Can be switched on and off. Computationally more expensive due to velocity field interpolation)
- Coupling via scatter-gather method
- Supports any number of diffusable species
- Supports general, customizable reaction terms (involving any number of interacting species and any number of kinetic coefficients)
- Dirichlet and Neumann boundary conditions
- Rectangular domains (static or expanding with constant speed)
- Chemotaxis mechanism 
- Helper functions to extract measurement data (readout position, precision-zone width, etc.)

## Output
A standard usage produces 3 output types:
- A series of ``.vtp`` frames containing the polygons to be visualized in paraview.
- A series of ``.vts`` frames containing the grid of the finite difference solver. Can be loaded into the same paraview session. 
- ``simulation.cfg`` saves all parameters (set in const.h) used for a simulation run for reproducibility. 

Additionally, depending on the experiment:
- A ``.off`` file for saving the polygon ensemble at a certain time point (usually at the end) to later load as the input/starting point of another simulation. 

## Folder Structure
`root`: makefile, input file (ensemble.off), README.md
`/include`: Contains all header files which make up the core of this software  
`/out`: Default output folder for files  

## Source files
In the following the most important components of the software are briefly explained. 

- ``const.h``: Contains all parameters and some settings (like enabling advection-dilution or chemotaxis). Treat this like a configuration-file. 

- ``domain.h``: Defines a very simple rectangular domain within which all calculations will happen. The boundary conditions of the solver will be applied at the boundaries of this domain. A constant growth factor can be set for each direction to turn the static domain into an expanding one.  

- ``ensemble.h``: The ensemble struct contains the core of the original PolyHoop software. It is responsible for the mechanical part of the simulation. It represents biological cells as polygons whose dynamics are governed by their potential energy.  

- ``geometry.h``: Contains the geometric structs used by the ensemble: Point, Vertex, Polygon. 

- ``grid.h``: Defines a simple matrix-like grid data structure used for the solver fields. Under the hood it is a simple nested std::vector. 

- ``solver.h``: The solver handles the reaction-diffusion part of the simulation. It uses central finite difference approximations and explicit euler time stepping. 
Additionally this file defines the available boundary conditions (Dirichlet or Neumann) which can be set individually for each side of the rectangular domain. By default all boundaries are treated as zero-flux.

- ``interpolator.h``: The interpolator takes care of moving data back and forth between the ensemble (polygons) and the solver (grid). It does this with a scatter-gather process similar to the particle-in-cell (PIC) method otherwise commonly used in plasma physics. Additionally, a velocity field is interpolated using the movement of the cells, more precisely their vertices.     

- ``utils.h``: Various quality-of-life functions.
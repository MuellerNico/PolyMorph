# PolyMorph

For the officially released version in [Computer Physics Communications](https://www.sciencedirect.com/science/article/pii/S0010465525000840), please see the branch [release-cpc](https://github.com/MuellerNico/PolyMorph/tree/release-cpc)

## Abstract

PolyMorph is an extension of [PolyHoop](https://www.sciencedirect.com/science/article/pii/S0010465524000511), a 2D mechanical tissue simulation software for soft particle dynamics developed by Vetter et al. PolyMorph adds a numerical solver for reaction-advection-diffusion equations which allows the modeling of chemical signaling. 

The simulation software is programmed in C++. It makes use of finite-difference approximations, explicit Euler time integration, scatter-gather interpolation between Eulerian and Lagrangian reference frames and shared memory parallelism with OpenMP.

Bachelor's thesis project of Nicolas Müller, ETH Zürich.

## Installation and execution
Use the provided makefile to compile the program. Then run the executable, specifying the number of cores available on your machine:
```shell
$ cd PolyMorph/
$ make
$ OMP_NUM_THREADS=4 ./polymorph
```
On a laptop with 4 cores @ 2.80GHz the default setup takes about 1 minute to run. PolyMorph does not use any external libraries. It only requires a compiler compatible with the C++11 standard and the OpenMP 3.1 specification for multi-threading support. 

Optionally, one can specify an output folder (default is `/out`):
```shell
$ OMP_NUM_THREADS=4 ./polymorph my/output/folder
```

## Directory structure
`root`: main.cpp, makefile, ensemble.off (input file), README.md    
`/include`: Contains all header files which make up the core of this software  
`/out`: Default output folder. Already contains the example output produced by the default setup when running the program out-of-the-box.  

## Source files

``main.cpp``: Simulation entry point and main time stepping loop. 

``param.h``: Contains all parameters and settings. Treat this like a configuration file. 

``domain.h``: Defines a simple rectangular domain within which all calculations happen. The boundary conditions of the solver will be applied at the boundaries of this domain.

``ensemble.h``: The ensemble struct contains the core of the original PolyHoop software. It is responsible for the mechanical part of the simulation. It represents biological cells as polygons whose dynamics are governed by their potential energy.  

``geometry.h``: Contains the geometric structs used by the ensemble: Point, Vertex, Polygon. 

``grid.h``: Defines a simple matrix-like grid data structure used for the FDM solver. Under the hood, a grid is simply a nested std::vector. 

``solver.h``: The solver handles the chemical reaction-diffusion part of the simulation. It uses central finite difference approximations and explicit euler time stepping. Additionally, this file defines the available boundary conditions (Dirichlet or Neumann) which can be set individually for each chemical species and each side of the rectangular domain. By default all boundaries are treated as zero-flux.

``interpolator.h``: The interpolator takes care of moving data back and forth between the ensemble (Lagrangian polygons) and the solver (Eulerian grid). It does this with a scatter-gather approach. Additionally, a velocity field is interpolated from the vertex velocities. 

``utils.h``: Various quality-of-life functions.

## Output
A standard usage produces 3 types of output files:
- A series of ``.vtp`` frames representing the polygon ensemble to be visualized in ParaView. 
- A series of ``.vts`` frames representing the grid of the finite difference solver to be visualized in the same ParaView session. 
- ``simulation.cfg`` saves all parameters set in param.h in text format for reproducibility. 
- Optionally, a ``.off`` file saving the polygon ensemble at a certain time point (e.g. at the end of a simulation run) to later load as the input/starting point of another simulation. 

## Usage
To customize the simulation for your application scenario, set the desired parameters in ``param.h`` and adjust other properties such as boundary conditions, reaction term, concentration effects on cell behavior, etc. inside ``main.cpp``. The program was designed with the intent, that the user only has to work with these two files. 

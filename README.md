# PolyMorph

Extension of PolyHoop for tissue morphogenesis coupled to chemical signaling

Copyright (c) 2024-2025  
Nicolas Müller, <nicolmueller@ethz.ch>  
Roman Vetter, <vetterro@ethz.ch>  
ETH Zurich

## Installation and execution
Use the provided makefile to compile the program. Then run the executable, specifying the number of cores to use, an optional input ensemble file (default: `ensemble.off`), and an optional output directory (default: `output`):

```shell
$ make
$ OMP_NUM_THREADS=4 ./polymorph ensemble.off output
```

PolyMorph does not use any external libraries. It only requires a compiler compatible with the C++11 standard and the OpenMP 3.1 specification for multi-threading support.

## Directory structure
`root`: main.cpp, makefile, ensemble.off, LICENSE, README.md

`./include`: Contains all header files which make up the core of this software

`./output`: Default output folder. Already contains the output produced by the default setup when running the program out-of-the-box.

## Files
``main.cpp``: Simulation entry point and main time stepping loop. Includes the definition of the functional relationships between cells and the chemical species.

``makefile``: Makefile containing compilation options.

``ensemble.off``: Initial tissue geometry in the Object File Format. The supplied file contains just a single circular cell with unit radius.

``LICENSE``: License file.

``README.md``: This readme file.

``include/param.h``: Contains all simulation parameters. Modify this in addition to ``main.cpp`` to customize the simulation setup.

``include/ensemble.h``: The ensemble contains the core of the original [PolyHoop software](https://doi.org/10.1016/j.cpc.2024.109128). It represents biological cells as interacting deformable polygons.

``include/geometry.h``: Defines several geometric objects such as Point, Vertex, Polygon.

``include/grid.h``: Defines a regular 2D grid data structure used for the finite difference solver.

``include/solver.h``: The finite difference solver handling the reaction-advecion-diffusion part of the simulation.

``include/interpolator.h``: The interpolator takes care of moving data back and forth between the ensemble (Lagrangian polygons) and the solver (Eulerian grid) with a scatter-gather approach. Additionally, the velocity field is interpolated from the polygon vertices.

``include/utils.h``: Various utility functions.

## Output
A standard usage produces the following output files:  

- A series of ``.vtp`` frames representing the polygon ensemble (the cellular tissue) for visualization in ParaView.
- A series of ``.vts`` frames representing the grid of the finite difference solver for visualization in ParaView.
- A text file named ``simulation.cfg`` listing all simulation parameters set in param.h for reproducibility.
- Optionally, a ``.off`` file saving the polygon ensemble at a certain time point (e.g. at the end of a simulation) to later load it as the input/starting point of another simulation.

## Usage
To customize the simulation for a specific application, set the desired parameters in ``param.h`` and adjust other properties such as boundary conditions, reaction term, concentration effects on cell behavior, etc. in ``main.cpp``. The program was designed such that only these two files need to be adjusted.

The default setup simulates a growing tissue with a single cell secreting a morphogen that the other cells read out to differentiate into a quiescent state. Additionally, the cells migrate up gradient of a second diffusing chemical. On a laptop with 4 cores at 2.80GHz, it takes about 1 minute to run.

# SPH
Set Partitioning Heuristic (SPH) based on CPLEX MILP solver and created starting from the CFT [1] Heuristic implemented by Accorsi Luca and Cavaliere Francesco (available soon).

[1] *Caprara, A., Fischetti, M., & Toth, P. (1999). A Heuristic Method for the Set Covering Problem. Operations Research, 47(5), 730–743. [doi:10.1287/opre.47.5.730](https://doi.org/10.1287/opre.47.5.730)*

It took part to DIMACS 12th Implementation Challenge on Vehicle Routing Problems in conjuction with: 
1. Improved version of Helsgaun's LKH-3.  
    **Team name:** _Vavavuma!_  
    **Team members:** F. Cavaliere, M. Fischetti, K. Helsgaun  
    **Algorithm name:** _LKHSP_  
    **Repository:** [https://github.com/c4v4/LKH-SP](https://github.com/c4v4/LKH-SP)) 

2. L. Accorsi and D. Vigo FILO algorithm.  
    **Team name:** _Route2Bo_  
    **Team members:** L. Accorsi, F. Cavaliere, D. Vigo  
    **Algorithm name:** FSP4D  
    **Repository:** [https://github.com/falcopt/f4d](https://github.com/falcopt/f4d)


## Build
1. Create a target binary target directory:

        [sph]$ mkdir build
        

2. Create Makefile using Cmake

        [sph]$ cd build
        [build]$ cmake ..
    
    The current build options are available:

        [build]$ cmake .. [-DCMAKE_BUILD_TYPE=[Release|Debug]] # Build type
        [build]$ cmake .. [-DENABLE_VERBOSE=[ON|OFF]]  # Switch on/off prints

3. Build the binary and the single-header with make

        [build]$ make -j
        
The binary will be placed inside the chosen directory.
The merged header file containing all the project will be placed in "`sph/one_header_only/SPH.hpp`".

## Use as library
To use the algorithm in other C++ project, just include generated header (`SPH.hpp`) in your project.

The SPHeuristic class can be used as high-level interface to the algorithm.

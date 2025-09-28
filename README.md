# SBP

This repository contains the implementation code of the paper "[Bayesian Adaptive Sparse Copula](https://www.economics.utoronto.ca/mburda/papers/sparsity.pdf)" by [Martin Burda](https://www.economics.utoronto.ca/mburda/) and [Artem Prokhorov](https://sites.google.com/site/artembprokhorov) (2025).

The code is written in [Modern Fortran](https://fortran-lang.org/), and it is portable across any compiler and platform that supports the F2008 standard. The code is designed with scalable [parallelization](https://wvuhpc.github.io/Modern-Fortran/20-Parallel-Programming/index.html). If [MPI](https://github.com/open-mpi/ompi) is available and enabled with a compilation flag (e.g. -D_MPI for mpifort) then MPI-specific code and communication libraries will be accessed by the compiler, resulting in faster run time on multi-core processors or distributed compting nodes. If MPI is not available, the slower serial version of the program will be compiled and can be executed on a single processor with a longer run time. 

Auxiliary [Rmd](https://rmarkdown.rstudio.com/index.html) files for data download, plotting of graphs and figures in the paper, are also provided, along with a [Makefile](https://fortran-lang.org/en/learn/building_programs/build_tools/#) with rules for code compilation and execution.

## Code Files
1. *SBP.f90* - the main program
2. *global_data.f90* - declaration of data structures
3. *arrays.f90* - routines for data allocation 
4. *timers.f90* - routines for wallclock, CPU, and MPI run time
5. *functions.f90* - implementation functions and routines
6. *marginals.f90* - routines for sampling of marginal densities
7. *MPI_communications.f90* - MPI instructions
8. *rootfindmod.f90* - routine for finding a root of a function

## External Dependencies
Additionally, the code also uses the following external dependencies:

1. Module `random` for generating random draws from the Beta and Dirichlet densities, available [here](https://www.netlib.org/random/random.f90) (save to file *random.f90* in the project directory).

2. The function `normal_01_cdf` for efficient evaluation of a standard Gaussian cdf, available [here](https://people.math.sc.edu/Burkardt/f_src/prob/prob.f90) (save in module `normals` in file *normals.f90* in the project directory). 

3. Module `sorts` containing sorting routines, available [here](
https://www.mjr19.org.uk/IT/sorts/sorts.f90) (save in file *sorts.f90* in the project directory).

4. Routine `newuoa` that solves unconstrained optimization problems without using derivatives, available [here](https://www.zhangzk.net/software.html) (save in file *newuoa.f90* in the project directory).

## Data Processing
1. *_Data download.Rmd* - code for data download from Yahoo Finance
2. *_Cdensity.Rmd* - plots copula densities used in the paper
3. *_Sparsity.Rmd* - plots of posterior sparsity over a grid
4. *_VaR_ES.Rmd* - plots of VaR and ES for a portfolio of assets

## Compilation
Compiler directives with flags for for [NVidia HPC SDK](https://developer.nvidia.com/hpc-sdk) are provided in *Makefile*. The compiler suite can be installed on a Windows machine using [WSL2](https://ubuntu.com/desktop/wsl) with [Ubuntu Desktop](https://ubuntu.com/desktop). 

Compile and run the code with:

<span style="font-family: monospace;">make run</span>

Object, module, and executable files are removed with:

<span style="font-family: monospace;">make clean</span>

# Simulations
The files used in the Monte Carlo simulation of the paper are included in the folder Simulations. We used R files (subfolder R_sim) for generating the simulated data, plotting figures, and obtaining summary statistics, and Fortran files (subfolder Fortran_sim) for running the code implementation on the simulated data. Bash shell files for compilation are included.

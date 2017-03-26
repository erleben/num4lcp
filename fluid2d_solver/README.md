# Fluid LCP

This code is supplementary code for the paper

A Fast Linear Complementarity Problem (LCP) Solver for Separating Fluid-Solid Wall Boundary Conditions by Andersen, Niebe and Erleben, VRIPHYS 2017

Here follows step by step tips on getting up and running

- Install CUDA toolkit (if not already installed, we used version 8.0)
- Install CMake (We used 3.6.2 but older versions should work too)
- Clone or download the git repository

In your local workspace, do as follows

- mkdir build
- cd build
- ccmake ../  (press confiugre 'c' 3 times and then generate 'g')
- make all
- cd ..
- cd bin
- ./fluid_solver default.cfg

You can copy the cfg-file or make changes to it as you please. Most of the settings are self-explanatory. The two most interesting ones are the, host and use_lcp, these determines whether the code will run on a GPU and whether one will use s standard Poisson solver or the LCP solver for separating fluid-solid wall boundary conditions.





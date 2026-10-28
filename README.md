# MiCRM_stability.jl

This package provides an implementation of the generalised microbial consumer resource model as detailed in the paper [Generalised modelling reveals stabilising mechanisms in cross-feeding microbial communities](https://arxiv.org/abs/2610.31154).

### Instalation
To install the package in a julia environment run the line
'''
add https://github.com/CleggTom/MiCRM_stability.jl.git
'''
in the package manager interface

### Worflow
The basic workflow of the package is as follows
1. User defines the structure of the uptake and stochiometric matricies themselves, or with one of the functions in 'community_generation.jl'
2. Generate a 'Parameter' object using 'get_structural_params' and by supplying a set of corresponding 'exponential_params'
3. Compute Jacobian matrix using the 'jacobian!' function
From there one can evaluate stability or any other properties of interest from the jaccobian matrix. 

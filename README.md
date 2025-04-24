# Iterative Reconstruction for Cosmological X-ray Tomography

Here we provide MATLAB implementations of iterative reconstruction methods 
for imaging cosmic strings using Cosmic Microwave Background (CMB) data.  
These examples use simulated data.

## Project Description

We develop numerical methods for inverting the light ray transform in 
relation to imaging cosmic strings in a two-dimensional universe over time
from  CMB measurements. We investigate various iterative reconstruction 
methods and analyze the artefact phenomena. 
 
## Installation 
### Software language

       MATLAB 9.14 (R2023a)
       For those without access to MATLAB, Octave provides an alternative platform.  
       Note that these codes have not been tested in Octave. 

### Requirements
The MainDrivers require the following package:

    "IR tools: A MATLAB Package of Iterative Regularization"
    by Silvia Gazzola, Per Christian Hansen and James G. Nagy
    https://github.com/jnagy1/IRtools.git

    "AIR Tools II: Algebraic Iterative Reconstruction Methods, Improved Implementation"
    by Per Christian Hansen and Jakob Sauer Jorgensen
    https://github.com/jakobsj/AIRToolsII

## How to Use
For Experiment 1 in [1], the main drivers is:
    
    Driver_Ex1.m          Compares reconstructions for the cone beam tomography 
                  reconstruction problem.  Iterative reconstruction methods are used.             

To run these codes, open a MATLAB command window, and type 
     
     >> Driver_Ex1 [press enter/return]
 
     
### Contributors
        Julianne Chung, 
        Department of Mathematics, Emory University
        
        Lucas Onisk, 
        Department of Mathematics, Emory University
        
        Yiran Wang, 
        Department of Mathematics, Emory University
	
## Licensing

If you use this codes, you *must* cite the original authors:

       [1] Chung, Onisk, and Wang. "Iterative Reconstruction Methods For 
        Cosmological X-ray Tomography", SIAM Journal on Imaging Sciences, 2025.


[MIT](LICENSE)

## Acknowledgement

This work was partially supported by the National Science Foundation under grants DMS-2341843, DMS-2038118, and DMS-2205266. Any opinions, findings, conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

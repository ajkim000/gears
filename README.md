# Hydrodynamic Spin-Coupling of 2D Rotors

**Alison Kim** \\
Courant Institute of Mathematical Sciences, NYU

![Simulation Video](media/animation.gif)

## Overview

Numerical simulation in MATLAB of the fluid-mediated coupling between two rigid, disk-shaped rotors immersed in a 2D viscous incompressible fluid. The left rotor ('active') is driven at a prescribed angular velocity. The emergent angular velocity of the right rotor ('passive') is measured to quantify fluid-mediated spin coupling.

The system is solved using the rigid penalty immersed boundary (pIB) method as well as an incompressible finite-difference FFT-based Navier–Stokes fluid solver on a periodic domain with semi-implicit timestepping.  

## How to run
1. Set key parameters (such as geometry, fluid properties, grid size) in initialize.m  
2. Run main.m

## Outputs
- Omega2_all — passive rotor angular velocity over time  
- Torque_all — hydrodynamic torque history  
- PNG frames of particle visualization  
- Periodic .mat simulation snapshots

## References
  Kim, Y. & Peskin, C.S. (2016).
    "A penalty immersed boundary method for a rigid body in fluid."
    Physics of Fluids, 28(3), 033603.
    
  Battista, N.A., Strickland, W.C., & Miller, L.A. (2017).
    "IB2d: a Python and MATLAB implementation of the immersed 
    boundary method." Bioinspiration & Biomimetics, 12(3), 036003.


# OceananigansLagrangianFilter.jl

*A package to compute Lagrangian temporal filters of Oceananigans.jl simulation output*

The Lagrangian filter can be run in two ways

- **Online**, integrated into your Oceananigans simulation. This avoids the need to save data at wave-resolving resolution, but the available filter shapes are not as desirable. 
- **Offline**, run after your Oceananigans simulation (or, feasibly, on any simulation output worked into the same format as Oceananigans native output) on saved data. Data should be at a temporal resolution that resolves the high frequency motions to be filtered. Velocities and the fields to be filtered need to be saved. The post-processing filter step runs similarly to an Oceananigans simulation, using the Oceananigans infrastructure to solve the filtering PDEs. 
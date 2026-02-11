---
title: 'OceananigansLagrangianFilter.jl: '
tags:
  - julia
  - ocean
  - Lagrangian filtering
  
authors:
  - name: Lois Baker
    orcid: 0000-0003-2678-3691
    corresponding: true
    affiliation: "1"
  
affiliations:
 - name: School of Mathematics and Maxwell Institute for Mathematical Sciences, University of Edinburgh, Edinburgh, UK
   index: 1


date: 10 February 2026
bibliography: paper.bib
---
# Summary
Numerical fluid dynamical simulations of geophysical fluids, such as the ocean and the atmosphere, typically exhibit complex dynamics including vortices and waves. The analysis of these simulations often requires the high frequency component (usually the waves) to be removed from the output to allow analysis of the lower frequency component, which we refer to as the 'mean flow'. Alternatvely, the faster component may be of interest, and must be separated from the mean flow for analysis purposes. Eulerian temporal filtering (or equivalently 'time averaging') is most commonly used for this purpose -- at a given spatial location, a time average (perhaps with some weight function) can be found either by post-processing simulation output 'offline', or more commonly at simulation time by constructing the average cumulatively 'online'. However, it has long been recognised that a Lagrangian time average, in which the property of interest is averaged in the frame of a (moving) fluid parcel, defines a more physically relevant mean flow [@andrewsExactTheoryNonlinear1978]. Furthermore, Lagrangian filtering, whereby waves are defined by an intrinisic (frame of flow) frequency criterion can most effectively decompose waves and mean flows [@bakerLagrangianFilteringWave2025]. Lagrangian filtering is usually performed using particle tracking methods, but its use has so far been limited by the extensive computational and memory requirements of the calculation. 

In this package, we implement methods for Lagrangian filtering based on the solution of partial differential equations (PDEs), building on recent theoretical developments [@kafiabadComputingLagrangianMeans2023;@bakerLagrangianFilteringWave2025;@minzEfficientLagrangianAveraging2025]. Although these methods are entirely generalisable to different numerical fluid dynamical solvers, here we implement them for use with the ``Oceananigans.jl`` ocean modelling framework. ``Oceananigans.jl`` is a popular and flexible Julia-based open source software package for finite volume simulations of the nonhydrostatic and hydrostatic Boussinesq equations on CPUs and GPUs [@Oceananigans]. ``OceananigansLagrangianFilter.jl`` provides a user-friendly way for ``Oceananigans.jl`` users to analyse their simulations and gain new insights using Lagrangian filtering, either by post-processing their existing simulation output ('offline') or integrating the functionality into their simulations at run-time ('online'). 


# Statement of Need
<!-- A section that clearly illustrates the research purpose of the software and places it in the context of related work. This should clearly state what problems the software is designed to solve, who the target audience is, and its relation to other work.-->

Lagrangian filtering consists of finding, for some simulated scalar field $f(x,t)$, the generalised Lagrangian mean $\overline{f}^{\mathrm{L}}(\symbf{x},t)$ defined by
$$\overline{f}^{\mathrm{L}}(\bar{\symbf{\varphi}}(\symbf{a},t),t) = \int_{-\infty}^\infty G(t-s) f(\symbf{\varphi}(\symbf{a},s),s)\,\mathrm{d} s \,,\label{eq:fbarL}$$
where $\symbf{\varphi}(\symbf{a},t)$ is the flow map, defining the position of a fluid particle with label $\symbf{a}$ at time $t$, $G(t)$ is some weight function (or equivalently the impulse response of the filter), and $\bar{\symbf{\varphi}}(\symbf{a},t)$ is the mean flow map, defined by 
$$\bar{\symbf{\varphi}}(\symbf{a},t) = \int_{-\infty}^\infty G(t-s) \symbf{\varphi}(\symbf{a},s)\mathrm{d} s \,.\label{eq:phibar}$$ 

For comparison, the Eulerian mean is defined as
$$\bar{f}(\symbf{x},t) = \int_{-\infty}^\infty G(t-s) f(\symbf{x},s)\,\mathrm{d} s \,.\label{eq:fE}$$
While $\overline{f}^{\mathrm{L}}$ is the true generalised Lagrangian mean as it is defined at the mean position of a trajectory, we also define a spatial rearrangement $f^*$ that is instead defined on the trajectory itself
$$f^*(\symbf{\varphi}(\symbf{a},t),t) = \int_{-\infty}^\infty G(t-s) f(\symbf{\varphi}(\symbf{a},s),s)\,\mathrm{d} s \,,\label{eq:fstar}$$
and a map
$$\symbf{\Xi}(\symbf{\varphi}(\symbf{a},t),t) = \bar{\symbf{\varphi}}(\symbf{a},t)\,,\label{eq:Xi}$$
that maps between $f^*$ and $\overline{f}^{\mathrm{L}}$, such that $\overline{f}^{\mathrm{L}}(\symbf{\Xi}(\symbf{x},t),t) = f^*(\symbf{x},t)$.

``OceananigansLagrangianFilter.jl`` solves directly for $f^*$ for any specified variables $f$ (including velocity components), which may be sufficient if $\symbf{\Xi}$ is close to the identity, or for certain wave decompositions [@bakerLagrangianFilteringWave2025]. Optionally, it also solves for $\symbf{\Xi}$ and afterward performs a spatial interpolation to recover $\overline{f}^{\mathrm{L}}$. The weight function $G(t)$ is constructed from sums of (complex) exponentials, which allows closed PDE systems for $f^*$ and $\symbf{\Xi}$ [@minzEfficientLagrangianAveraging2025] (see \autoref{fig:filtershape}).

Lagrangian filtering is superior to Eulerian filtering in two primary situations. The first is when waves are Doppler-shifted by the mean flow, such that they can be low frequency in the rest frame of the simulation. The canonical example of this is internal lee waves, which are generated by steady stratified flow over topography. \autoref{fig:leewave}a shows the horizontal velocity field $u$ generated by a steady flow $U = 0.1$ ms$^{-1}$ with stratification $N = 0.001$ s$^{-1}$ over a Gaussian bump, which was run in 2D using the `HydrostaticFreeSurfaceModel` of `Oceananigans.jl`. 

The waves are phase-locked to the topography such that their frequency in the rest frame is zero. When an Eulerian average is taken, the waves therefore remain in the mean component as shown in \autoref{fig:leewave}b. However, internal gravity waves are subject to a criterion $\omega > f$ on their intrinsic frequency $\omega$, where $f$ is the Coriolis parameter, which means that they can be filtered from the mean flow using a Lagrangian filter that removes high frequencies. \autoref{fig:leewave}c shows $u^*$, demonstrating that the wave perturbations have been removed - note the order of magnitude difference in the colorbar range between the top and bottom rows. The effect of the wave displacement remains in $u^*$, as visualised by contours of the Lagrangian filtered buoyancy $b^*$ in \autoref{fig:leewave}c. However, when $u^*$ and $b^*$ are remapped to $\overline{u}^{\mathrm{L}}$ and $\overline{b}^{\mathrm{L}}$ respectively, the waves are entirely removed (\autoref{fig:leewave}d). 

![Here's a caption. \label{fig:IO}](figures/offline_IO_lagrangian_filtering.png)

The second use-case is when the wave (or other high frequency oscillation) is large compared to the scales of the mean flow -- for example, the displacement of a submesoscale eddy field by an inertial oscillation, or tidal displacements of coastal waters. Figures \ref{fig:IO}a,e show an idealised and non-dimensional `Oceananigans.jl` simulation using the `ShallowWaterModel`, in which a Gaussian blob of tracer $T$ is advected in a circular path by a spatially uniform inertial oscillation. The Eulerian average (Figures \ref{fig:IO}b,f) leads to blurring of the tracer gradients due to its high-frequency displacement. The Lagrangian average $T^*$ (Figures \ref{fig:IO}c,g) is in this case identical to the raw field $T$ (Figures \ref{fig:IO}a,e) since the tracer is conserved along a trajectory, but when remapped to $\overline{T}^{\mathrm{L}}$ (Figures \ref{fig:IO}d,h) the effect of the displacement is removed and the tracer maintains its original shape, as in the absence of the inertial oscillation. 

![Here's a caption. \label{fig:leewave}](figures/offline_lee_wave_lagrangian_filtering.png)

Explain here briefly the benefit of being able to make this separation (cite other filtering papers) - identification of internal waves in simulations. 

# State of the field
<!-- A description of how this software compares to other commonly-used packages in the research area. If related tools exist, provide a clear “build vs. contribute” justification explaining your unique scholarly contribution and why existing alternatives are insufficient. -->
- Generally particle tracking, both online and offline. Online methods often suffer from clustering and interpolation errors, offline from memory requirements. In particular, it is generally necessary to track particles along an entire trajectory to find the Lagrangian mean at one timestep.
- Also explain the online/offline nature and why the flexibility is good.  
- This method allows calculation of the Lagrangian mean/and or high frequency component at every timestep. Talk about relation to Shakespeare 2021 particle tracking. 


# Software design
<!-- An explanation of the trade-offs you weighed, the design/architecture you chose, and why it matters for your research application. This should demonstrate meaningful design thinking beyond a superficial code structure description.-->
- Several options for implementation of Lagrangian mean via PDEs - why I chose this one and offline and online options. 
- Offline and online (include figure of filter shapes) different uses in different contexts. Keep as much of Oceananigans infrastructure as possible, building on the existing constructors.
- Refer to equations? To demonstrate what is actually found by different methods, cite Baker 25 for f*. It is a choice to find f* first (no interpolation) 

![Caption for filter shape. \label{fig:filtershape}](figures/filter_responses.png)

Some more text here
# Research impact statement
<!-- Evidence of realized impact (publications, external use, integrations) or credible near-term significance (benchmarks, reproducible materials, community-readiness signals). The evidence should be compelling and specific, not aspirational. -->
- Makes available for implementation in other fluid dynamical solvers
- Analytical methods have been published, but this implementation allows the methods to be easily accessible to the community. 
- Useful for calculating Lagrangian means for saving as transport velocities

# AI usage disclosure
<!-- Transparent disclosure of any use of generative AI in the software creation, documentation, or paper authoring. If no AI tools were used, state this explicitly. If AI tools were used, describe how they were used and how the quality and correctness of AI-generated content was verified.-->

# Acknowledgements
Thank NFFDy. We want to thank the [Climate Modeling Alliance](https://clima.caltech.edu) team and ``Oceananigans.jl`` contributors for their fantastic project. We also thank ``JuliaGPU`` contributors for ``KernelAbstractions.jl`` and ``CUDA.jl`` which make this code possible. 


# References
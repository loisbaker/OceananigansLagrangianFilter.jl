# Online Filtering

This page describes the Lagrangian filtering equations and their implementation in the 'online' configuration of *OceananigansLagrangianFilter.jl*. They are very similar to the forward pass of the 'offline' configuration (described in [Offline Filtering](@ref "Offline Filtering")).

## Equations
We directly compute
```math
\begin{equation}\label{fstardef}
    f^*(\vb*{\varphi}(\vb*{a},t),t) = \int_{-\infty}^t G(t-s)f(\vb*{\varphi}(\vb*{a},s),s) \, \mathrm{d} s\,,
\end{equation}
```
and optionally compute
```math
\begin{equation}\label{Xidefonline}
\vb*{\Xi}(\vb*{\varphi}(\vb*{a},t),t) = \int_{-\infty}^t \alpha e^{-\alpha(t-s)}\vb*{\varphi}(\vb*{a},s)\,\mathrm{d} s\,,
\end{equation}
```
so that the generalised Lagrangian mean ``\bar{f}^{\mathrm{L}}`` (see definition in [Lagrangian averaging](@ref "Lagrangian averaging")) can be recovered by a post-processing interpolation step using
```math
\begin{equation}
\bar{f}^{\mathrm{L}}(\vb*{\Xi}(\vb*{x},t),t) = f^*(\vb*{x},t)\,.
\end{equation}
```

We consider filter kernels composed of sums of ``N`` exponentials, where ``N`` is even.
```math
\begin{equation}\label{onlinekernel}
    G(t) = \begin{cases}
    \sum_{n=1}^{N/2} e^{-c_n |t|} \left( a_n \cos(d_n |t|) + b_n \sin(d_n |t|) \right)\,, \hspace{1cm} t > 0\,,\\
    0\,,\hspace{1cm} t \leq 0\,.
    \end{cases}
\end{equation}
```
We require ``G`` to be normalised such that
```math
\begin{equation}
    \int_{-\infty}^\infty G(s) \, ds = 1\,,
\end{equation}
```
or equivalently, that ``\hat{G}(0) = 1``. This requires that
```math
\begin{equation}
    \sum_{n=1}^{N/2} \frac{a_nc_n + b_n d_n}{c_n^2 + d_n^2} = 1\,.
\end{equation}
```
This normalisation is only strictly required when we define a map that computes the trajectory mean position (`map_to_mean = true` in [`OnlineFilterConfig`](@ref "OnlineFilterConfig")) but we keep the requirement for now.

We define a set of ``N`` weight functions. For ``k = 1,...,N/2`` we have
```math
\begin{align}
    G_{Ck}(t) &= \begin{cases}
    e^{-c_k t}\cos d_k t\,,\hspace{1cm} t > 0\,, \\
    0\,,\hspace{1cm} t \leq 0\,,
    \end{cases}
    G_{Sk}(t) &= \begin{cases}
    e^{-c_k t}\sin d_k t\,, \hspace{1cm} t > 0\,,\\
    0\,,\hspace{1cm} t \leq 0\,.
    \end{cases}
\end{align}
```
For ``t>0``, we have
```math
\begin{align} \label{forward_G_derivs}
    G_{Ck}'(t) &= - c_kG_{Ck}(t) - d_k G_{Sk}(t)\,, \\
    G_{Sk}'(t) &= - c_kG_{Sk}(t) + d_k G_{Ck}(t)\,. \\
\end{align}
```

We then define a corresponding set of ``N`` filtered scalars
```math
\begin{align}\label{forwardgdef}
    g_{Ck}(\vb*{\varphi}(\vb*{a},t),t) &= \int_{-\infty}^t G_{Ck}(t-s)f(\vb*{\varphi}(\vb*{a},s),s)\, \mathrm{d} s\,,\\
    g_{Sk}(\vb*{\varphi}(\vb*{a},t),t) &= \int_{-\infty}^t G_{Sk}(t-s)f(\vb*{\varphi}(\vb*{a},s),s)\, \mathrm{d} s\,,\\
\end{align}
```
so that 
```math
\begin{equation}\label{reconstitutefstar}
    f^*(\vb*{x},t) = \sum_{n=1}^{N/2} a_n g_{Cn}(\vb*{x},t) + b_n g_{Sn}(\vb*{x},t)\,.
\end{equation}
```
We derive PDEs for the ``g_{Ck}`` and ``g_{Sk}`` by taking the time derivative of \eqref{forwardgdef}:
```math
\begin{align}
\frac{\partial g_{Ck}}{\partial t} + \vb*{u} \cdot \nabla g_{Ck} &= f - c_k g_{Ck} - d_k g_{Sk} \label{gCeqn}\\
\frac{\partial g_{Sk}}{\partial t} + \vb*{u} \cdot \nabla g_{Sk} &=  - c_k g_{Sk} + d_k g_{Ck} \label{gSeqn}
\end{align}
```
This system of equations can be solved with initial conditions ``g_{Ck}(\vb*{x},0) = g_{Sk}(\vb*{x},0)=0`` (TODO add more on ICs, spin-up)

If we want to find ``\bar{f}^{\mathrm{L}}``, we define map functions with which to interpolate ``f^*`` to ``\bar{f}^{\mathrm{L}}`` after the simulation. We define 
```math
\begin{align}
    \vb*{\Xi}_{Ck}(\vb*{\varphi}(\vb*{a},t),t) &= \int_{-\infty}^t G_{Ck}(t-s)\vb*{\varphi}(\vb*{a},s) \mathrm{d} s\\
    \vb*{\Xi}_{Sk}(\vb*{\varphi}(\vb*{a},t),t) &= \int_{-\infty}^t G_{Sk}(t-s)\vb*{\varphi}(\vb*{a},s) \mathrm{d} s\,,
\end{align}
```
so that
```math
\begin{equation}
\vb*{\Xi}(\vb*{x},t) = \sum_{n=1}^{N/2} a_n \vb*{\Xi}_{Ck}(\vb*{x},t) + b_n \vb*{\Xi}_{Sk}(\vb*{x},t)
\end{equation}
```
and
```math
\begin{align}
\frac{\partial \vb*{\Xi}_{Ck}}{\partial t} + \vb*{u} \cdot \nabla \vb*{\Xi}_{Ck} &= \vb*{x} - c_k \vb*{\Xi}_{Ck} - d_k \vb*{\Xi}_{Sk} \\
\frac{\partial \vb*{\Xi}_{Sk}}{\partial t} + \vb*{u} \cdot \nabla \vb*{\Xi}_{Sk} &=  - c_k \vb*{\Xi}_{Sk} + d_k \vb*{\Xi}_{Ck}\,.
\end{align}
```
To define perturbation equations, we set:
```math
\begin{align}
    \vb*{\Xi}_{Ck} &= \vb*{\xi}_{Ck} + \frac{c_k}{c_k^2 + d_k^2}\vb*{x} \\
    \vb*{\Xi}_{Sk} &= \vb*{\xi}_{Sk} + \frac{d_k}{c_k^2 + d_k^2}\vb*{x}\,,
\end{align}
```
where the coefficients of ``\vb*{x}`` are needed because each of the filters ``G_{Ck}`` and ``G_{Sk}`` are not individually normalised over the interval ``[-\infty,0]``.
The perturbation map equations are then given by
```math
\begin{align}
\frac{\partial \vb*{\xi}_{Ck}}{\partial t} + \vb*{u} \cdot \nabla \vb*{\xi}_{Ck} &= -\frac{c_k}{c_k^2 + d_k^2}\vb*{u} - c_k \vb*{\xi}_{Ck} - d_k \vb*{\xi}_{Sk} \label{xiCeqn}\\
\frac{\partial \vb*{\xi}_{Sk}}{\partial t} + \vb*{u} \cdot \nabla \vb*{\xi}_{Sk} &= -\frac{d_k}{c_k^2 + d_k^2}\vb*{u} - c_k \vb*{\xi}_{Sk} + d_k \vb*{\xi}_{Ck}\,,\label{xiSeqn}
\end{align}
```
and solved with initial conditions ``\vb*{\xi}_{Ck}(\vb*{x},0) = \vb*{\xi}_{Sk}(\vb*{x},0)=0``.

## Implementation

The online Lagrangian filter equations are solved at the same time as the governing equations of the simulation to be filtered. This means that equations \eqref{gCeqn} and \eqref{gSeqn}, and optionally \eqref{xiCeqn} and \eqref{xiSeqn} if remapping is required, need to be passed to the `model` that is being used. 

OceananigansLagrangianFilter provides helper functions to define these extra fields and their forcings, and an example is given in [geostrophic_adjustment.jl](https://github.com/loisbaker/OceananigansLagrangianFilter.jl/blob/main/examples/geostrophic_adjustment.jl). 

Other helper functions are provided to initialise the filtered variables, define output fields that compute ``f^*`` from \eqref{reconstitutefstar}, regrid ``f^*`` to ``\bar{f}^{\mathrm{L}}``, compute the Eulerian filter for comparison, compute a shifted time variable, and output a final netcdf file. 

We summarise here the protocol for implementing online Lagrangian filtering

- Install OceananigansLagrangianFilter in your environment (see [Quick Start](@ref "Quick Start")).

- Load OceananigansLagrangianFilter and its utility functions at the top of your script. This automatically loads Oceananigans.

```julia
using OceananigansLagrangianFilter
using OceananigansLagrangianFilter.Utils
```

- Setup your parameters, grid, tracers, and forcing as normal

- Define your `filter_config` - an [`OnlineFilterConfig`](@ref "OnlineFilterConfig"). This takes as arguments:
    - `grid`: the grid you have already defined
    - `output_filename`: a filename to save filtered output to
    - `var_names_to_filter`: a tuple of strings defining the names of variables to filter. These can be any of your `tracers`. Velocities to filter don't need to be listed here (see below).
    - `velocity_names`: the velocity names that you want to use to compute Lagrangian trajectories. These are also the velocities that will be filtered if `compute_mean_velocities = true`.
    - `N` and `freq_c`: Can be provided together to give a Butterworth filter of order ``N`` with cutoff frequency `freq_c`. 
    - `filter_params`: a named tuple of coefficients `a1`, `b1`, `c1`, `d1`, `a2`, `b2`, `c2`, `d2`, etc defining a filter kernel as in \eqref{onlinekernel}
    - `map_to_mean`: A Bool determining whether to compute the maps ``\vb*{\xi}_{Ck}`` and ``\vb*{\xi}_{Sk}`` and solve \eqref{xiCeqn} and \eqref{xiSeqn}.
    - `compute_mean_velocities`: A Bool determining whether to compute and output the mean velocities. They are computed from the maps ``\vb*{\xi}_{Ck}`` and ``\vb*{\xi}_{Sk}``, so if `map_to_mean=false` and `compute_mean_velocities=true` the maps will still be computed. 

```julia
filter_config = OnlineFilterConfig( grid = grid,
                                    output_filename = filename_stem * ".jld2",
                                    var_names_to_filter = ("b","T"), 
                                    velocity_names = ("u","w"),
                                    N = 2,
                                    freq_c = f/2)
```

- Create the filtered variables ``g_{Ck}``, ``g_{Sk}``, ``\vb*{\xi}_{Ck}`` and ``\vb*{\xi}_{Sk}`` using the function [`create_filtered_vars`](@ref "create_filtered_vars"), which only needs the `filter_config`. These filtered variables will be added to the model as tracers.

```julia
filtered_vars = create_filtered_vars(filter_config)
```

- Create forcing for these filtered variables using the function [`create_forcing`](@ref "create_forcing"). This implements the right-hand-sides of equations \eqref{gCeqn}, \eqref{gSeqn}, \eqref{xiCeqn}, and \eqref{xiSeqn}. Merge this `filter_forcing` with your existing forcing. 

```julia
filter_forcing = create_forcing(filtered_vars, filter_config)
forcing = merge(forcing, filter_forcing);
```

- If you're using a closure, you'll need to tell the filtered scalars not to use a closure (although they could be given a closure if necessary for stability - this has proved unecessary so far and is more accurate). A helper function [`zero_closure_for_filtered_vars`](@ref "zero_closure_for_filtered_vars") to set the diffusivity to zero for each of the filtered variables is provided.

```julia
zero_filtered_var_closure = zero_closure_for_filtered_vars(filter_config)
```

- Define the closure for your model variables using the `zero_filtered_var_closure`

```julia
horizontal_closure = HorizontalScalarDiffusivity(ν=1e-6, κ=merge((T=1e-6, b= 1e-6),zero_filtered_var_closure) )
vertical_closure = VerticalScalarDiffusivity(ν=1e-6 , κ=merge((T=1e-6, b= 1e-6),zero_filtered_var_closure) )
closure = (horizontal_closure, vertical_closure)
```

- Define the model with the combined (original and filtered) `tracers`, `forcing` and `closure`. This should work with both `NonHydrostaticModel` and `HydrostaticFreeSurfaceModel`.

- Initialise your model variables as normal

- Initialise the filtered variables (this uses the previously initialised model variables to give a better initialisation, so needs to be performed after the previous step)

```julia
initialise_filtered_vars_from_model(model, filter_config)
```

- Define the simulation, any callbacks, `conjure_time_step_wizard`, etc as normal

- Use the [`create_output_fields`](@ref "create_output_fields") helper function to define the output fields (this defines the outputs ``f^*`` and ``\vb*{\Xi}`` so that all of the intermediate filter variables are not output by default, though they could be examined as for any other tracer)

```julia
outputs = create_output_fields(model, filter_config)
```

- Add in any other output fields, including the original fields if desired (not included by default in the online filter).

```julia
outputs["b"] = model.tracers.b
outputs["T"] = model.tracers.T
outputs["u"] = model.velocities.u
outputs["v"] = model.velocities.v
outputs["w"] = model.velocities.w
```

- Define a `JLD2Writer` for the `outputs`. The filename should be the same as `filter_config.output_filename`. For the post-processing functions provided in the next few steps, this does need to be a `JLD2Writer` rather than a `NetCDFWriter`, but a helper function is provided to output a final `.nc` file if needed. 

- Run the simulation
    
```julia
run!(simulation)
```
- Optionally regrid to mean position using [`regrid_to_mean_position!`](@ref "regrid_to_mean_position!"). This adds a new field to the output data file.

```julia
if filter_config.map_to_mean
    regrid_to_mean_position!(filter_config)
end
```

- Optionally compute the Eulerian filter (with the same `filter_params`) using [`compute_Eulerian_filter!`](@ref "compute_Eulerian_filter!").
    
```julia
compute_Eulerian_filter!(filter_config)
```

- Optionally compute a shifted time coordinate to give a more appropriate reference time for the average. The new reference time is calculated as the weighted mean of time: ``t_{shift} = \int_{-\infty}^t G(t-s)s\, \mathrm{d} s``
    
```julia
compute_time_shift!(filter_config)
```

- Optionally output a final NetCDF file using [`jld2_to_netcdf`](@ref "jld2_to_netcdf"). This helper function should work for any `.jld2` Oceananigans output. 

```julia
jld2_to_netcdf(filename_stem * ".jld2", filename_stem * ".nc")
```
# Implementation of online filtering

The online Lagrangian filter equations, which find Lagrangian filtered tracer(s) ``f^*`` (see [Lagrangian averaging](@ref "Lagrangian averaging") for a definition) are solved at the same time as the governing equations of your simulation, in which the original tracer ``f`` is being found. Filtered velocities can also be found. This means that equations for the filtered tracers, and optionally the maps (see [Online Lagrangian filtering equations](@ref "Online Lagrangian filtering equations")) need to be passed to your `model`. 

OceananigansLagrangianFilter provides helper functions to define these extra fields and their forcings, and an example is given in [geostrophic_adjustment.jl](https://github.com/loisbaker/OceananigansLagrangianFilter.jl/blob/main/examples/online_filter_geostrophic_adjustment.jl). 

Other helper functions are provided to initialise the filtered variables, define output fields that compute ``f^*``, regrid ``f^*`` to ``\bar{f}^{\mathrm{L}}`` (see [Lagrangian averaging](@ref "Lagrangian averaging")), compute the Eulerian filter for comparison, compute a shifted time variable, and output a final NetCDF file. 

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
    - `filter_params`: a named tuple of coefficients `a1`, `b1`, `c1`, `d1`, `a2`, `b2`, `c2`, `d2`, etc defining a filter kernel (see [Choosing online filters](@ref "Choosing online filters")).
    - `map_to_mean`: A Bool determining whether to compute the maps ``\vb*{\xi}_{Ck}`` and ``\vb*{\xi}_{Sk}`` and solve their equations (see [Online Lagrangian filtering equations](@ref "Online Lagrangian filtering equations")).
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

- Create forcing for these filtered variables using the function [`create_forcing`](@ref "create_forcing"). This implements the right-hand-sides of the tracer and map equations (see [Online Lagrangian filtering equations](@ref "Online Lagrangian filtering equations")). Merge this `filter_forcing` with your existing forcing. 

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
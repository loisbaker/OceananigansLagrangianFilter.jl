# Online filtering implementation

The online Lagrangian filter equations, which find Lagrangian filtered tracer(s) ``f^*`` (see [Lagrangian averaging](@ref "Lagrangian averaging") for a definition) are solved at the same time as the governing equations of your simulation, in which the original tracer ``f`` is being found. Filtered velocities can also be found. This means that equations for the filtered tracers, and optionally the maps (see [Online Lagrangian filtering equations](@ref "Online Lagrangian filtering equations")) need to be passed to your `model`. 

OceananigansLagrangianFilter provides helper functions to define these extra fields and their forcings, and an example is given at [Geostrophic adjustment online](@ref "Geostrophic adjustment with online Lagrangian filtering").

Other helper functions are provided to initialise the filtered variables, define output fields that compute ``f^*``, regrid ``f^*`` to ``\bar{f}^{\mathrm{L}}`` (see [Lagrangian averaging](@ref "Lagrangian averaging")), compute the Eulerian filter for comparison, compute a shifted time variable, and output a final NetCDF file. 

We summarise here the protocol for implementing online Lagrangian filtering

- Install OceananigansLagrangianFilter in your environment (see [Quick Start](@ref "Quick Start")).

- Load OceananigansLagrangianFilter and its utility functions at the top of your script. This automatically loads Oceananigans.

```julia
using OceananigansLagrangianFilter
using OceananigansLagrangianFilter.Utils
```

- Setup your parameters, grid, tracers, and forcing as normal

- Define your `filter_config` - an [`OnlineFilterConfig`](@ref). This takes as arguments:
    - `grid`: the grid you have already defined
    - `output_filename`: a filename to save filtered output to
    - `var_names_to_filter`: a tuple of strings defining the names of variables to filter. These can be any of your `tracers`. Velocities to filter don't need to be listed here (see below).
    - `velocity_names`: the velocity names that you want to use to compute Lagrangian trajectories. These are also the velocities that will be filtered if `compute_mean_velocities = true`.
    - `N` and `freq_c`: Can be provided together to give a Butterworth filter of order ``N`` with cutoff frequency `freq_c`. 
    - `filter_params`: a named tuple of coefficients `a1`, `b1`, `c1`, `d1`, `a2`, `b2`, `c2`, `d2`, etc defining a filter kernel (see [Choosing online filters](@ref "Choosing online filters")).
    - `compute_maps`: A Bool determining whether to compute and output the maps ``\vb*{\xi}_{Ck}`` and ``\vb*{\xi}_{Sk}`` by solving their equations (see [Online Lagrangian filtering equations](@ref "Online Lagrangian filtering equations")). The maps are the displacements from the mean position.
    - `regrid_to_mean`: A Bool determining whether the filtered fields should be interpolated to the mean position after the simulation (see below). This requires the maps, so if `regrid_to_mean=true` the maps will be computed.
    - `compute_mean_velocities`: A Bool determining whether to compute and output the mean velocities. They are computed from the maps ``\vb*{\xi}_{Ck}`` and ``\vb*{\xi}_{Sk}``, so if `compute_maps=false` and `compute_mean_velocities=true` the maps will still be solved for (but not output). 
    - `lagrangian`: `true` (default) for the Lagrangian filter, or `false` for the Eulerian filter (see below).
    - `relaxation_mask_func` and `relaxation_params`: optionally, relaxation of the filtered variables towards their
      local equilibrium where the mask is non-zero, e.g. near open boundaries, as for the offline filter (see
      [Boundary relaxation](@ref "Boundary relaxation")).
    - `frequency_scale_func`, `frequency_scale_params` and `frequency_scale_dependencies`: optionally, a frequency scale that varies the cutoff frequency in space and time, as for the offline filter (see [Varying cutoff frequency](@ref "Varying cutoff frequency")). Online, the dependencies can be any fields of the model.

```julia
filter_config = OnlineFilterConfig( grid = grid,
                                    output_filename = filename_stem * ".jld2",
                                    var_names_to_filter = ("b","T"), 
                                    velocity_names = ("u","w"),
                                    N = 2,
                                    freq_c = f/2)
```

- Create the filtered variables ``g_{Ck}``, ``g_{Sk}``, ``\vb*{\xi}_{Ck}`` and ``\vb*{\xi}_{Sk}`` using the function [`create_filtered_vars`](@ref), which only needs the `filter_config`. These filtered variables will be added to the model as tracers.

```julia
filtered_vars = create_filtered_vars(filter_config)
```

- Create forcing for these filtered variables using the function [`create_forcing`](@ref). This implements the right-hand-sides of the tracer and map equations (see [Online Lagrangian filtering equations](@ref "Online Lagrangian filtering equations")). Merge this `filter_forcing` with your existing forcing. 

```julia
filter_forcing = create_forcing(filtered_vars, filter_config)
forcing = merge(forcing, filter_forcing);
```

- If you're using a closure, you'll need to tell the filtered scalars not to use a closure (although they could be given a closure if necessary for stability - this has proved unecessary so far and is more accurate). A helper function [`zero_closure_for_filtered_vars`](@ref) to set the diffusivity to zero for each of the filtered variables is provided.

```julia
zero_filtered_var_closure = zero_closure_for_filtered_vars(filter_config)
```

- Define the closure for your model variables using the `zero_filtered_var_closure`

```julia
horizontal_closure = HorizontalScalarDiffusivity(ν=1e-6, κ=merge((T=1e-6, b= 1e-6),zero_filtered_var_closure) )
vertical_closure = VerticalScalarDiffusivity(ν=1e-6 , κ=merge((T=1e-6, b= 1e-6),zero_filtered_var_closure) )
closure = (horizontal_closure, vertical_closure)
```

- For the Lagrangian filter (the default), the filtered variables are advected like any other tracer. For the Eulerian filter (`lagrangian = false`), they must not be advected: set the model's `tracer_advection` with the helper function [`filtered_vars_advection`](@ref), which gives no advection (`nothing`) for each of the filtered variables. Tracers missing from a `NamedTuple` `tracer_advection` are advected with `Centered()`, not the model's `advection`, so list all the other tracers too. If the advection of the filtered variables doesn't match `lagrangian`, [`initialise_filtered_vars_from_model`](@ref) and [`create_output_fields`](@ref) give a warning.

```julia
tracer_advection = merge((; b = WENO(), T = WENO()), filtered_vars_advection(filter_config))
```

- Define the model with the combined (original and filtered) `tracers`, `forcing` and `closure` (and `tracer_advection` for the Eulerian filter). This should work with both `NonHydrostaticModel` and `HydrostaticFreeSurfaceModel`. Any fields named in `frequency_scale_dependencies` must be fields of the model, e.g. added to its `auxiliary_fields` (which Oceananigans doesn't update, so keep them up to date, e.g. with a callback).

- Initialise your model variables as normal

- Initialise the filtered variables (this uses the previously initialised model variables to give a better initialisation, so needs to be performed after the previous step)

```julia
initialise_filtered_vars_from_model(model, filter_config)
```

- Define the simulation, any callbacks, `conjure_time_step_wizard`, etc as normal

- Use the [`create_output_fields`](@ref) helper function to define the output fields (this defines the outputs ``f^*`` and ``\vb*{\Xi}`` so that all of the intermediate filter variables are not output by default, though they could be examined as for any other tracer)

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
- Optionally regrid to mean position using [`regrid_to_mean_position!`](@ref). This adds a new field to the output data file.

```julia
if filter_config.regrid_to_mean
    regrid_to_mean_position!(filter_config)
end
```

- Optionally compute the Eulerian filter (with the same `filter_params`) using [`compute_Eulerian_filter!`](@ref). With a frequency scale, this is only available if it depends on the coordinates alone (no `frequency_scale_dependencies`); otherwise, see `lagrangian = false` above.
    
```julia
compute_Eulerian_filter!(filter_config)
```

- Optionally compute a shifted time coordinate to give a more appropriate reference time for the average. The new reference time is calculated as the weighted mean of time: ``t_{shift} = \int_{-\infty}^t G(t-s)s\, \mathrm{d} s``. This isn't available with a frequency scale.
    
```julia
compute_time_shift!(filter_config)
```

- Optionally output a final NetCDF file using [`jld2_to_netcdf`](@ref). This helper function should work for any `.jld2` Oceananigans output. 

```julia
jld2_to_netcdf(filename_stem * ".jld2", filename_stem * ".nc")
```
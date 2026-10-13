# Offline filtering implementation

The offline Lagrangian filter equations, which find Lagrangian filtered tracer(s) ``f^*`` (see [Lagrangian averaging](@ref "Lagrangian averaging") for a definition) are solved after the original Oceananigans simulation (or, feasibly, using any simulation output worked into the same format as Oceananigans native output) on saved data. Data should be at a temporal resolution that resolves the high frequency motions to be filtered. Velocities and the tracer(s) ``f`` to be filtered need to be provided. The post-processing filter step runs similarly to an Oceananigans simulation, using the Oceananigans infrastructure to solve the filtering PDEs. 

The offline filter uses mostly the same functions as the online filter to define filtered fields and their forcings, but in this case most of the process is 'under the hood', as the user only needs to provide the simulation data and specify the configuration. An example is given in [`offline_filter_geostrophic_adjustment.jl`](https://github.com/loisbaker/OceananigansLagrangianFilter.jl/blob/main/examples/offline_filter_geostrophic_adjustment.jl), and more detail is given in [how it works](@ref "How it works").

## Input data

The original data can be a JLD2 or NetCDF file of Oceananigans output, given by `original_data_filename`. It must contain the velocities (`velocity_names`) and the variables to be filtered (`var_names_to_filter`), saved often enough to resolve the motions being filtered.

- **JLD2** files written by Oceananigans' `JLD2Writer` can be used directly.
- **NetCDF** files should be written by Oceananigans' `NetCDFWriter` (or have the same metadata), since the grid and the location of each variable are read from the file's attributes. The time coordinate must be called `time`.

The data is read frame by frame directly from the original file, so memory use doesn't depend on the length of the dataset, and no intermediate copy of the data is written (see [How it works](@ref "How it works")).

## Example

A short example of how to implement offline filtering on a GPU is given below:

```julia
using OceananigansLagrangianFilter
using Oceananigans.Units
using CUDA

# Define the filter configuration
filter_config = OfflineFilterConfig(original_data_filename = "my_simulation.jld2", # Where the original simulation output is
                                    output_filename = "my_filtered_simulation.jld2" # Where to save the filtered output
                                    var_names_to_filter = ("T", "b"), # Which variables to filter
                                    velocity_names = ("u","v"), # Velocities to use for Lagrangian filtering
                                    architecture = GPU(), # CPU() or GPU()
                                    Δt = 20minutes, # Time step of filtering simulation
                                    T_out = 1hour, # How often to output filtered data
                                    N = 2, # Order of Butterworth filter
                                    freq_c = 1e-4/2, # Cut-off frequency of Butterworth filter
                                    output_netcdf = false, # Whether to output filtered data to a netcdf file in addition to .jld2
                                    delete_intermediate_files = true, # Delete the individual output of the forward and backward passes
                                    compute_mean_velocities = true, # Whether to compute the mean velocities
                                    compute_Eulerian_filter = true) # Whether to compute the Eulerian filter for comparison

# Run the offline filter
run_offline_Lagrangian_filter(filter_config)

# The filtered data is now saved to `my_filtered_simulation.jld2`
```

## Frequency scale

The cutoff frequency can vary in space, and in time, through a frequency scale ``M`` that multiplies the filter's
frequencies locally (see [Varying cutoff frequency](@ref)). `frequency_scale_func` takes the non-flat coordinates and
`frequency_scale_params`:

```julia
frequency_scale(x, z, p) = 1 + 0.5 * sin(2π * x / p.L)

config = OfflineFilterConfig(original_data_filename = "my_simulation.jld2",
                             var_names_to_filter = ("T",),
                             velocity_names = ("u", "w"),
                             N = 2, freq_c = 1e-4,
                             frequency_scale_func = frequency_scale,
                             frequency_scale_params = (; L = 10kilometers))
```

``M`` can also depend on other fields in the data, which are read in with `extra_auxiliary_var_names` and named in
`frequency_scale_dependencies`. Their values at the cell centre are passed after the coordinates:

```julia
frequency_scale(x, z, zeta, p) = 1 + abs(zeta) / p.f

config = OfflineFilterConfig(...,
                             extra_auxiliary_var_names = ("zeta",),
                             frequency_scale_func = frequency_scale,
                             frequency_scale_dependencies = (:zeta,),
                             frequency_scale_params = (; f = 1e-4))
```

- ``M`` must be finite and positive. This is checked at every cell centre when the config is built, but with
  dependencies only for the initial fields, when the filtered variables are initialised. Make sure it stays
  positive throughout the run, e.g. by using `abs` or `exp` of the fields rather than the fields themselves.
- `frequency_scale_func` is called inside the filter's kernels, so it must be GPU-compatible to run on a GPU.
- Fields in `frequency_scale_dependencies` must be at cell centres, apart from the velocities, which are interpolated
  to cell centres.
- ``M`` can't depend on the filtered fields of the same run (see [Varying cutoff frequency](@ref)): to use a
  filtered field, filter it in a first run and read it in with `extra_auxiliary_var_names` in a second.
- Choose `Δt` for the largest local cutoff during the run, and leave enough time for spin-up at the smallest one.
- The time shift is not defined, and the Eulerian filter in post-processing is not computed. For the Eulerian
  counterpart, run the filter again with `advection = nothing`.

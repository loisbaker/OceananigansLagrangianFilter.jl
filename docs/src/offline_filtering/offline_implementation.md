# Offline filtering implementation

The offline Lagrangian filter equations, which find Lagrangian filtered tracer(s) ``f^*`` (see [Lagrangian averaging](@ref "Lagrangian averaging") for a definition) are solved after the original Oceananigans simulation (or, feasibly, using any simulation output worked into the same format as Oceananigans native output) on saved data. Data should be at a temporal resolution that resolves the high frequency motions to be filtered. Velocities and the tracer(s) ``f`` to be filtered need to be provided. The post-processing filter step runs similarly to an Oceananigans simulation, using the Oceananigans infrastructure to solve the filtering PDEs. 

The offline filter uses mostly the same functions as the online filter to define filtered fields and their forcings, but in this case most of the process is 'under the hood', as the user only needs to provide the simulation data and specify the configuration. An example is given in [`offline_filter_geostrophic_adjustment.jl`](https://github.com/loisbaker/OceananigansLagrangianFilter.jl/blob/main/examples/offline_filter_geostrophic_adjustment.jl), and more detail is given in [how it works](@ref "How it works").

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
## Spatially varying cutoff frequency

The offline filter accepts a stationary, positive `cutoff_mask` that sets the
filter-clock rate `dτ/dt = cutoff_mask` along trajectories. A constant mask
value `m` gives the usual local cutoff `m * freq_c`; a varying mask adapts as
particles move. See [Spatially varying cutoff equations](@ref) for the kernel.

```julia
input = FieldTimeSeries("my_simulation.jld2", "T")
cutoff_mask = Field{Center, Nothing, Center}(input.grid)
set!(cutoff_mask, (x, z) -> 1 + 0.5sin(2π * x / 10_000))

config = OfflineFilterConfig(original_data_filename = "my_simulation.jld2",
                             var_names_to_filter = ("T",),
                             velocity_names = ("u", "v", "w"),
                             grid = input.grid, N = 2, freq_c = 1e-4,
                             cutoff_mask = cutoff_mask)
```

The mask must be finite and strictly positive. `Center` or `Nothing` locations
allow invariant dimensions. A scalar mask uses the existing uniform-cutoff path.
For boundary relaxation, set `boundary_relaxation`, `relax_timescale`,
`mask_func`, and `mask_params` in the same config. `mask_func` selects the
relaxation region; `cutoff_mask` sets the local filter and map-equilibrium target.
The optional Eulerian post-processing filter does not yet accept spatial masks.

Choose `Δt` for the largest local cutoff and leave enough time for endpoint
transients at the smallest one. Keep the mask fixed during a run.

# Offline filtering implementation

The offline Lagrangian filter equations, which find Lagrangian filtered tracer(s) ``f^*`` (see [Lagrangian averaging](@ref "Lagrangian averaging") for a definition) are solved after the original Oceananigans simulation (or, feasibly, using any simulation output worked into the same format as Oceananigans native output) on saved data. Data should be at a temporal resolution that resolves the high frequency motions to be filtered. Velocities and the tracer(s) ``f`` to be filtered need to be provided. The post-processing filter step runs similarly to an Oceananigans simulation, using the Oceananigans infrastructure to solve the filtering PDEs. 

The offline filter uses mostly the same functions as the online filter to define filtered fields and their forcings, but in this case most of the process is 'under the hood', as the user only needs to provide the simulation data and specify the configuration. An example is given in [offline_filter_geostrophic_adjustment.jl](https://github.com/loisbaker/OceananigansLagrangianFilter.jl/blob/main/examples/offline_filter_geostrophic_adjustment.jl), and more detail on how it works is given in [How it works](@ref "How it works").

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
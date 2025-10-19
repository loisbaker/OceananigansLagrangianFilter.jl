# 🚀 Quick Start

### Offline Filtering

Offline filtering (whereby the data is processed after simulation time) allows for better filter shapes, since for a given reference time, data from the past and the future is available. The filters implemented here have real frequency response, and therefore have linear phase shift. If the exact properties of the filter shape are important, then offline filtering is preferable. 

Here is a simple example of how to filter a pre-existing dataset.

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
You can find an example of a simple simulation of geostrophic adjustment followed by offline filtering in [`offline_filter_geostrophic_adjustment.jl`](https://github.com/loisbaker/OceananigansLagrangianFilter.jl/blob/main/examples/offline_filtering_geostrophic_adjustment.jl).

### Online Filtering

For online filtering, you would integrate the filter directly into your `Oceananigans.jl` setup, using the helper functions provided. This is explained full in [Online filtering implementation](@ref "Online filtering implementation"), and an example is given in [`online_filtering_geostrophic_adjustment.jl`](https://github.com/loisbaker/OceananigansLagrangianFilter.jl/blob/main/examples/online_filtering_geostrophic_adjustment.jl). The filtered values are then computed as your simulation runs, avoiding the need to save data at high frequency. 

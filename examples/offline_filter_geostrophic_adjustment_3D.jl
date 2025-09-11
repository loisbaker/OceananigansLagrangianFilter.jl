using OceananigansLagrangianFilter
using Oceananigans.Units
using CairoMakie
using CUDA

# Set up the filter configuration
filter_config = OfflineFilterConfig(original_data_filename="geostrophic_adjustment_3D.jld2", # Where the original simulation output is
                                    output_filename = "geostrophic_adjustment_filtered_3D.jld2", # Where to save the filtered output
                                    var_names_to_filter = ("T", "b"), # Which variables to filter
                                    velocity_names = ("u","v","w"), # Velocities to use for Lagrangian filtering
                                    architecture = GPU(), # CPU() or GPU()
                                    Δt = 20minutes, # Time step of filtering simulation
                                    T_out = 1hour, # How often to output filtered data
                                    N = 2, # Order of Butterworth filter
                                    freq_c = 1e-4/2, # Cut-off frequency of Butterworth filter
                                    output_netcdf = true, # Whether to output filtered data to a netcdf file in addition to .jld2
                                    delete_intermediate_files = true, # Delete the individual output of the forward and backward passes
                                    compute_Eulerian_filter = true) # Whether to compute the Eulerian filter for comparison


# Run the offline Lagrangian filter
run_offline_Lagrangian_filter(filter_config)


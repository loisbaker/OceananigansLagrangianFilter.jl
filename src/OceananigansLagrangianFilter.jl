module OceananigansLagrangianFilter

using Reexport
@reexport using Oceananigans  

# Define a supertype for all configuration objects.
abstract type AbstractConfig end

# Define submodules
include("./Utils/Utils.jl")
include("./OfflineLagrangianFilter/OfflineLagrangianFilter.jl")
include("./OnlineLagrangianFilter/OnlineLagrangianFilter.jl")

using .Utils
using .OfflineLagrangianFilter
using .OnlineLagrangianFilter

# User will have access to these functions directly from the main module
export OfflineFilterConfig, OnlineFilterConfig, run_offline_Lagrangian_filter

# User will also have direct access to utility functions
export create_input_data_on_disk, load_data, set_offline_BW2_filter_params, set_online_BW_filter_params
export create_original_vars, create_filtered_vars, create_forcing, create_output_fields
export update_input_data!, initialise_filtered_vars_from_data, initialise_filtered_vars_from_model, zero_closure_for_filtered_vars
export sum_forward_backward_contributions!, regrid_to_mean_position!
export jld2_to_netcdf, get_weight_function, get_frequency_response, compute_Eulerian_filter!, compute_time_shift!

end # module
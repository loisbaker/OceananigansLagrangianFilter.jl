module Utils

using OceananigansLagrangianFilter: AbstractConfig
using JLD2
using JLD2: Group
using Oceananigans
using Oceananigans.Fields: Center
using Oceananigans.Units: Time

# Need to figure out how to get OfflineLagrangianFilter
include("lagrangian_filter_utils.jl")
include("post_processing_utils.jl")

export create_input_data_on_disk, load_data, set_offline_BW2_filter_params, set_online_BW_filter_params,
export create_original_vars, create_filtered_vars, create_forcing, create_output_fields
export update_input_data!, initialise_filtered_vars_from_data, initialise_filtered_vars_from_model, zero_closure_for_filtered_vars
export sum_forward_backward_contributions!, regrid_to_mean_position!
export jld2_to_netcdf, get_weight_function, get_frequency_response, compute_Eulerian_filter!, compute_time_shift!

end
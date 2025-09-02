module OfflineLagrangianFilter

using Reexport
@reexport using Oceananigans   

using DocStringExtensions

using KernelAbstractions: @index, @kernel

using Oceananigans.Utils
using Oceananigans.Grids
using Oceananigans.Solvers

using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: reconstruct_global_grid, Distributed
using Oceananigans.Grids: XYRegularRG, XZRegularRG, YZRegularRG, XYZRegularRG
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Utils: sum_of_velocities

using JLD2
using JLD2: Group
using Oceananigans.Fields: Center
using Oceananigans.Units: Time

import Oceananigans: fields, prognostic_fields 
import Oceananigans.Advection: cell_advection_timescale
import Oceananigans.OutputWriters: default_included_properties

# Need to look at which of these really need to be exported
export create_input_data_on_disk, load_data, set_BW_filter_params, create_original_vars, create_filtered_vars, create_forcing, create_output_fields, update_input_data!, sum_forward_backward_contributions!, regrid_to_mean_position!,jld2_to_netcdf
export c_div_U
export default_included_properties

include("lagrangian_filter.jl")
include("compute_lagrangian_filter_buffer_tendencies.jl")
include("compute_lagrangian_filter_tendencies.jl")
include("../LagrangianFilterUtils/lagrangian_filter_utils.jl")
include("../LagrangianFilterUtils/post_processing_utils.jl")
include("lagrangian_filter_tendency_kernel_functions.jl")
include("lagrangian_filtering_advection_operators.jl")
include("set_lagrangian_filter.jl")
include("show_lagrangian_filter.jl")
include("update_lagrangian_filter_state.jl")


#####
##### Update some Oceananigans internal methods for our new OfflineLagrangianFilter model.
#####

function cell_advection_timescale(model::OfflineLagrangianFilter)
    grid = model.grid
    velocities = model.velocities
    return cell_advection_timescale(grid, velocities)
end

"""
    fields(model::LagrangianFilter)

Return a flattened `NamedTuple` of the fields in `model.velocities`, `model.tracers`, and any
auxiliary fields for a `LagrangianFilter` model.
"""
fields(model::OfflineLagrangianFilter) = merge(model.velocities,
                                           model.tracers,
                                           model.auxiliary_fields)

"""
    prognostic_fields(model::LagrangianFilter)

Return a flattened `NamedTuple` of the prognostic fields associated with `LagrangianFilter`.
"""
prognostic_fields(model::OfflineLagrangianFilter) = merge(model.velocities, model.tracers)

#####
##### Define a config structure
#####

"""
    OfflineFilterConfig(;
        time_step = 1.0,
        filter_order = 2,
        output_filename = "offline_filtered_data.jld2",
        verbose = true,
    )

A configuration object for `apply_offline_filter`.
"""
Base.@kwdef struct OfflineFilterConfig
    original_data_filename::String = joinpath(@__DIR__, "my_simulation.jld2")
    T_start::Union{Real, Nothing} = nothing
    T_end::Union{Real, Nothing} = nothing
    arch::AbstractArchitecture = nothing
    T_out::Real = 0.1
    N::Int = 2
    freq_c::Real = 1.0
    var_names_to_filter::Tuple{Vararg{String}} = (,)
    velocity_names::Tuple{Vararg{String}} = (,)
    filter_params::Union{NamedTuple, Nothing} = nothing
    Δt::Real = 1e-3
    map_to_mean::Bool = true
    forward_output_filename::String = "forward_LF.jld2"
    backward_output_filename::String = "backward_LF.jld2"
    combined_output_filename::String = "combined_LF.jld2"

    verbose::Bool = true
    # Add any other configurable parameters here
end

#####
##### Run the filtering operation
#####
function do_offline_Lagrangian_filter(config)

    # Copy and manipulate data on disk to have correct order and time shift
    T = create_input_data_on_disk(config.original_data_filename, config.var_names_to_filter, config.velocity_names; direction = "forward", T_start = config.T_start, T_end = config.T_end)

  # Do the rest 


export do_offline_Lagrangian_filter
end # module OfflineLagrangianFilter


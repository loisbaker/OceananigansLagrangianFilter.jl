module OceananigansLagrangianFilter

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


import Oceananigans: fields, prognostic_fields 
import Oceananigans.Advection: cell_advection_timescale
import Oceananigans.OutputWriters: default_included_properties

using JLD2
using JLD2: Group
using Oceananigans
using Oceananigans.Fields: Center
using Oceananigans.BoundaryConditions: PeriodicBoundaryCondition
using Oceananigans.Units: Time

# We need a python import to use the LinearNDInterpolator for regridding
using PythonCall
const scipy_interpolate = PythonCall.pynew()
function __init__()
    PythonCall.pycopy!(scipy_interpolate, pyimport("scipy.interpolate"))
end

using ProgressBars
export set_data_on_disk!, load_data, set_BW_filter_params, create_original_vars, create_filtered_vars, create_forcing, create_output_fields, update_input_data!, sum_forward_backward_contributions!, regrid_to_mean_position!,jld2_to_netcdf
export c_div_U
export default_included_properties

include("LagrangianFilter/lagrangian_filter.jl")
include("LagrangianFilter/compute_lagrangian_filter_buffer_tendencies.jl")
include("LagrangianFilter/compute_lagrangian_filter_tendencies.jl")
include("LagrangianFilter/lagrangian_filter_utils.jl")
include("LagrangianFilter/lagrangian_filter_tendency_kernel_functions.jl")
include("LagrangianFilter/lagrangian_filtering_advection_operators.jl")
include("LagrangianFilter/set_lagrangian_filter.jl")
include("LagrangianFilter/show_lagrangian_filter.jl")
include("LagrangianFilter/update_lagrangian_filter_state.jl")



#####
##### AbstractModel interface
#####

function cell_advection_timescale(model::LagrangianFilter)
    grid = model.grid
    velocities = model.velocities
    return cell_advection_timescale(grid, velocities)
end

"""
    fields(model::LagrangianFilter)

Return a flattened `NamedTuple` of the fields in `model.velocities`, `model.tracers`, and any
auxiliary fields for a `LagrangianFilter` model.
"""
fields(model::LagrangianFilter) = merge(model.velocities,
                                           model.tracers,
                                           model.auxiliary_fields)

"""
    prognostic_fields(model::LagrangianFilter)

Return a flattened `NamedTuple` of the prognostic fields associated with `LagrangianFilter`.
"""
prognostic_fields(model::LagrangianFilter) = merge(model.velocities, model.tracers)


# This is a now-defunct part of OutputWriters/jld2_writer.jl that is needed in the earlier version of Oceananigans being temporarily used
default_included_properties(model::LagrangianFilter) = [:grid]

export LagrangianFilter
end # module OceananigansLagrangianFilter

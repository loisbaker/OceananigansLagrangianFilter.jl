module Utils

using OceananigansLagrangianFilter: AbstractConfig, AbstractOfflineConfig, AbstractOnlineConfig
using JLD2
using JLD2: Group
using Oceananigans
using Oceananigans.Fields: Center, location, indices
using Oceananigans.Units: Time
using Oceananigans.OutputReaders: InMemory, load_serialized_grid
using Oceananigans.OutputWriters: serializeproperty!
using Oceananigans: AbstractModel, boundary_conditions
using Oceananigans.BoundaryConditions: PeriodicBoundaryCondition, FieldBoundaryConditions, fill_halo_regions!
using Oceananigans.Grids: AbstractGrid, topology, Flat, on_architecture, λnodes, φnodes, xspacings, yspacings, zspacings
using Oceananigans.ImmersedBoundaries: immersed_cell
using KernelAbstractions: @kernel, @index
using DataStructures: OrderedDict
using NCDatasets
using Oceananigans.Forcings: DiscreteForcing
using Oceananigans.Grids: node
using Oceananigans.Operators: interpolation_operator, assumed_field_location

using ..DataIO: BufferedDataReader, create_buffered_reader, advance_buffer!, interpolate_to_model!, update_input_data!, stored_times, shared_halo_regions
using ..DataIO: create_data_source, cpu_frame_templates, read_frame!

include("filter_params.jl")
include("filter_variables.jl")
include("filter_forcing.jl")
include("filter_outputs.jl")
include("initialisation.jl")
include("file_utils.jl")
include("post_processing.jl")
include("regrid_to_mean.jl")

export create_buffered_reader, set_offline_BW2_filter_params, set_online_BW_filter_params, set_offline_spectrum_filter_params, set_online_spectrum_filter_params
export create_original_vars, create_filtered_vars, create_forcing, create_output_fields
export update_input_data!, initialise_filtered_vars_from_data, initialise_filtered_vars_from_model, zero_closure_for_filtered_vars, filtered_vars_advection
export sum_forward_backward_contributions!, regrid_to_mean_position!, change_sign_of_map_variables!
export jld2_to_netcdf, get_weight_function, get_frequency_response, compute_Eulerian_filter!, compute_time_shift!

end
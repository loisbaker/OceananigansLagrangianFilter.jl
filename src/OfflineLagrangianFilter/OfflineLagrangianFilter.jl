module OfflineLagrangianFilter

using DocStringExtensions
using KernelAbstractions: @index, @kernel

using Oceananigans
using Oceananigans.Utils
using Oceananigans.Grids
using Oceananigans.Solvers
using JLD2

using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: reconstruct_global_grid, Distributed
using Oceananigans.Grids: XYRegularRG, XZRegularRG, YZRegularRG, XYZRegularRG
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Utils: sum_of_velocities
using Oceananigans.OutputReaders: AbstractInMemoryBackend
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Architectures

import Oceananigans: fields, prognostic_fields 
import Oceananigans.Advection: cell_advection_timescale, AbstractAdvectionScheme
using ..OceananigansLagrangianFilter: AbstractConfig

export OfflineFilterConfig, run_offline_Lagrangian_filter, LagrangianFilter

using ..Utils

include("run_offline_lagrangian_filter.jl")
include("lagrangian_filter.jl")
include("compute_lagrangian_filter_buffer_tendencies.jl")
include("compute_lagrangian_filter_tendencies.jl")
include("lagrangian_filter_tendency_kernel_functions.jl")
include("lagrangian_filtering_advection_operators.jl")
include("set_lagrangian_filter.jl")
include("show_lagrangian_filter.jl")
include("update_lagrangian_filter_state.jl")

#####
##### Update some Oceananigans internal methods for our new LagrangianFilter model.
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

Return a flattened `NamedTuple` of the prognostic fields associated with `LagrangianFilter`. velocities are not really prognostic, but we include them here for convenience.
"""
prognostic_fields(model::LagrangianFilter) = merge(model.velocities, model.tracers)

#####
##### Define a config structure
#####

"""
    OfflineFilterConfig(;
        
    )

A configuration object for `apply_offline_filter`.
"""
struct OfflineFilterConfig <: AbstractConfig
    original_data_filename::String
    var_names_to_filter::Tuple{Vararg{String}}
    velocity_names::Tuple{Vararg{String}}
    T_start::Real
    T_end::Real
    T::Real
    architecture::AbstractArchitecture 
    T_out::Real 
    filter_params::NamedTuple
    Δt::Real
    backend::AbstractInMemoryBackend
    map_to_mean::Bool
    forward_output_filename::String
    backward_output_filename::String
    output_filename::String
    npad::Int
    delete_intermediate_files::Bool
    compute_Eulerian_filter::Bool
    output_netcdf::Bool
    output_original_data::Bool
    advection::AbstractAdvectionScheme
    grid::AbstractGrid
    filter_mode::String

end

"""
    OfflineFilterConfig(; original_data_filename::String,
                        var_names_to_filter::Tuple{Vararg{String}},
                        velocity_names::Tuple{Vararg{String}},
                        T_start::Union{Real,Nothing} = nothing,
                        T_end::Union{Real,Nothing} = nothing,
                        T::Union{Real,Nothing} = nothing,
                        architecture::AbstractArchitecture = CPU(),
                        T_out::Union{Real,Nothing} = nothing,
                        N::Union{Int, Nothing} = nothing,
                        freq_c::Union{Int, Nothing} = nothing,
                        filter_params::Union{NamedTuple, Nothing} = nothing,
                        Δt::Union{Real,Nothing} = nothing,
                        backend::AbstractInMemoryBackend = InMemory(4),
                        map_to_mean::Bool = true,
                        forward_output_filename::String = "forward_output.jld2",
                        backward_output_filename::String = "backward_output.jld2",
                        output_filename::String = "filtered_output.jld2",
                        npad::Int = 5,
                        delete_intermediate_files::Bool = true,
                        compute_Eulerian_filter::Bool = false,
                        output_netcdf::Bool = false,
                        output_original_data::Bool = true
                        advection::AbstractAdvectionScheme = WENO(),
                        grid::Union{AbstractGrid, Nothing} = nothing,
                        filter_mode::String = "offline")

Constructs a configuration object for offline Lagrangian filtering of Oceananigans data.
This function validates the input data file, time specifications, and filter parameters
before creating the `OfflineFilterConfig` object.

Keyword arguments
=================

  - `original_data_filename`: (required) The path to the JLD2 file containing the original Oceananigans output data.
  - `var_names_to_filter`: (required) A `Tuple` of `String`s specifying the names of the tracer variables to be filtered.
  - `velocity_names`: (required) A `Tuple` of `String`s specifying the names of the velocity fields in the data file to be used for advection.
  - `T_start`: Start time for the filter. Must be within the time range of the data. If not given, defaults to either T_end - T (if they are given), or the start time of the original data.
  - `T_end`: End time for the filter. Must be within the time range of the data. If not given, defaults to either T_start + T (after T_start given or computed, if T is given), or the end time of the original data.
  - `T`: Duration of the filtering. If not given, defaults to T_end - T_start (after T_start and T_end are given or computed).
  - `architecture`: The architecture (CPU or GPU) to be used for the filtering computation. Default: `CPU()`.
  - `T_out`: The output time step for the filtered data. If `nothing`, it defaults to the time step of the original data.
  - `N`, `freq_c`: Parameters for a Butterworth squared filter. `2^N` is the order of the filter, and `freq_c` is the cutoff frequency. These are used to automatically generate `filter_params` if not provided. Must be specified together if `filter_params` is not given.
  - `filter_params`: A `NamedTuple` containing the coefficients for a custom filter. Only filter_params OR `N` and `freq_c` should be given.
  - `Δt`: The time step for the internal Lagrangian filter simulation. If `nothing`, it defaults to `T_out / 10`, but this may not be appropriate.
  - `backend`: The backend for loading `FieldTimeSeries` data. See `Oceananigans.Fields.FieldTimeSeries`. Default: `InMemory(4)`.
  - `map_to_mean`: A `Bool` indicating whether to map filtered data to the mean position (i.e. calculate generalised Lagrangian mean). Default: `true`.
  - `forward_output_filename`: The filename for the output of the forward filter pass. Default: `"forward_output.jld2"`.
  - `backward_output_filename`: The filename for the output of the backward filter pass. Default: `"backward_output.jld2"`.
  - `output_filename`: The filename for the final combined and mapped output. Default: `"filtered_output.jld2"`.
  - `npad`: The number of cells to pad the interpolation to mean position, used when there are periodic boundary conditions. Default: `5`.
  - `delete_intermediate_files`: A `Bool` indicating whether to delete `forward_output.jld2` and `backward_output.jld2` after the final combined file is created. Default: `true`.
  - `compute_Eulerian_filter`: A `Bool` indicating whether to also compute an Eulerian-mean-based filter for comparison. Default: `false`.
  - `output_netcdf`: A `Bool` indicating whether to also convert the final JLD2 output file to a NetCDF file. Default: `false`.
    - `output_original_data`: A `Bool` indicating whether to include the original data in the final output file for comparison. Default: `true`.
  - `advection`: The advection scheme to use for the Lagrangian filter simulation. Default: `WENO()`. Using lower-order schemes may be a source of error.
  - `grid`: The grid for the simulation. If `nothing`, the grid is inferred from the `original_data_filename` (preferred option)
  - `filter_mode`: A `String` indicating whether to run the filter in "offline" or "online" mode. Default: "offline". TODO use multiple dispatch for this instead.
"""
function OfflineFilterConfig(; original_data_filename::String,
                            var_names_to_filter::Tuple{Vararg{String}},
                            velocity_names::Tuple{Vararg{String}},
                            T_start::Union{Real,Nothing} = nothing,
                            T_end::Union{Real,Nothing} = nothing,
                            T::Union{Real,Nothing} = nothing,
                            architecture::AbstractArchitecture = CPU(),
                            T_out::Union{Real,Nothing} = nothing,
                            N::Union{Int, Nothing} = nothing,
                            freq_c::Union{Real, Nothing} = nothing,
                            filter_params::Union{NamedTuple, Nothing} = nothing,
                            Δt::Union{Real,Nothing} = nothing,
                            backend::AbstractInMemoryBackend = InMemory(4),
                            map_to_mean::Bool = true,
                            forward_output_filename::String = "forward_output.jld2",
                            backward_output_filename::String = "backward_output.jld2",
                            output_filename::String = "filtered_output.jld2",
                            npad::Int = 5,
                            delete_intermediate_files::Bool = true,
                            compute_Eulerian_filter::Bool = false,
                            output_netcdf::Bool = false,
                            output_original_data::Bool = true,
                            advection::AbstractAdvectionScheme = WENO(),
                            grid::Union{AbstractGrid, Nothing} = nothing)

    # Check that the original file exists 
    if !isfile(original_data_filename)
        error("Source file not found: $original_data_filename")
    end

    # Open up the file for some checks
    jldopen(original_data_filename,"r") do original_file

        # Check if velocities and tracers given are in the original_file
        for var in (var_names_to_filter..., velocity_names...)
            if !haskey(original_file["timeseries"], var)
                error("Variable '$var' not found in original data file.")
            end
        end

        # Check for consistent T_start, T_end, T
        iterations = parse.(Int, keys(original_file["timeseries/t"]))
        times = [original_file["timeseries/t/$iter"] for iter in iterations]

        # If all are specified, make sure they're consistent and if not, throw an error
        if isnothing(T_start) + isnothing(T_end) + isnothing(T) == 0
            if T_start + T != T_end
                error("Inconsistent time specifications: T_start + T != T_end")
            elseif T-start > T_end
                error("Inconsistent time specifications: T_start > T_end")
            end

        # If 0,1, or 2 are specified, calculate the others
        elseif isnothing(T_start)
            if !isnothing(T_end) && !isnothing(T)
                T_start = T_end - T
            elseif !isnothing(T_end) && isnothing(T)
                T_start = times[1]
                T = T_end - T_start
            elseif !isnothing(T) && isnothing(T_end)
                T_start = times[1]
                T_end = T_start + T
            else # isnothing(T) && isnothing(T_end)
                T_start = times[1]
                T_end = times[end]
                T = T_end - T_start
            end
        elseif isnothing(T_end)
            if isnothing(T)
                T_end = times[end]
                T = T_end - T_start
            else # !isnothing(T)
                T_end = T_start + T
            end
        else # !isnothing(T)
            T = T_end - T_start
            if T < 0
                error("Inconsistent time specifications: T_start > T_end")
            end
        end

        # Check that T_start and T_end are found in the original_file
        if T_start < times[1] || T_start > times[end]
            error("T_start=$T_start is outside the range of the original data: [$times[1], $times[end]].")
        end

        if T_end < times[1] || T_end > times[end]
            error("T_end=$T_end is outside the range of the original data: [$times[1], $times[end]].")
        end

        @info "Filter interval will be from T_start=$T_start to T_end=$T_end, duration T=$T"

        # Now check T_out, set if necessary to same as input
        if isnothing(T_out)
            T_out = times[2] - times[1]
            @info "T_out not set. Setting T_out = $T_out"
        end

        # If Δt not set, set to T_out/10
        if isnothing(Δt)
            Δt = T_out / 10
            @info "Δt (filter simulation timestep) not set. Setting Δt = $Δt, but be careful"
        end
    end

    # Make sure we have some filter parameters
    if !isnothing(filter_params) && (!isnothing(N) || !isnothing(freq_c))
        error("Specify either filter_params or N and freq_c, not both.")
    elseif isnothing(filter_params) && (isnothing(N) || isnothing(freq_c))
        error("Must specify either filter_params or both N and freq_c.")
    elseif isnothing(filter_params) && !isnothing(N) && !isnothing(freq_c)
        filter_params = set_offline_BW_filter_params(;N,freq_c)
        @info "Setting filter parameters to use Butterworth squared, order $(2^N), cutoff frequency $freq_c"
    else # !isnothing(filter_params) && isnothing(N) && isnothing(freq_c)
        # User has specified filter_params directly, but we should check it has the right fields
        if haskey(filter_params, :N_coeffs)
            if N_coeffs == 0.5 # Single exponential special case
                if !all(haskey(filter_params, :a1) , haskey(filter_params, :c1))
                    error("For N_coeffs=0.5, filter_params must have fields :N_coeffs, :a1, and :c1")
                end
            elseif floor(filter_params.N_coeffs) != filter_params.N_coeffs
                error("N_coeffs must be a positive integer or 0.5")
            else
                if !all(haskey(filter_params, Symbol(coeff,i)) for coeff in ["a","b","c","d"] for i in 1:filter_params.N_coeffs)
                    error("For N_coeffs>0.5, filter_params must have fields :N_coeffs, :a1, :a2, ..., :b1, :b2, ..., :c1, :c2, ..., :d1, :d2, ...")
                end
            
            end
        else # N_coeffs isn't provided, but we might be able to infer it
            if floor(length(filter_params)/4) == length(filter_params)/4
                filter_params = merge(filter_params, (N_coeffs = length(filter_params)/4,))
                # But we still have to check that the right entries are there:
                if !all(haskey(filter_params, Symbol(coeff,i)) for coeff in ["a","b","c","d"] for i in 1:filter_params.N_coeffs)
                    error("filter_params must have fields :N_coeffs, :a1, :a2, ..., :b1, :b2, ..., :c1, :c2, ..., :d1, :d2, ...")
                end
            elseif length(filter_params) == 2
                filter_params = merge(filter_params, (N_coeffs = 0.5,))
                if !all(haskey(filter_params, :a1) , haskey(filter_params, :c1))
                    error("For a filter with two coefficients, filter_params must have fields :N_coeffs, :a1, and :c1")
                end
            else
                error("filter_params must have either 2 entries (for single exponential) or 4*N entries (for Butterworth squared of order 2^N)")
            
            end

        end
    end

    # Check normalisation of filter coefficients
    if filter_params.N_coeffs == 0.5
        if filter_params.a1/2 != filter_params.c1
            @warn "Filter coefficients are not normalised: a1=$(filter_params.a1) != c1=$(filter_params.c1)"
        end
    else
        a_coeffs = [filter_params[Symbol("a",i)] for i in 1:filter_params.N_coeffs]
        b_coeffs = [filter_params[Symbol("b",i)] for i in 1:filter_params.N_coeffs]
        c_coeffs = [filter_params[Symbol("c",i)] for i in 1:filter_params.N_coeffs] 
        d_coeffs = [filter_params[Symbol("d",i)] for i in 1:filter_params.N_coeffs]
        if sum((a_coeffs.*c_coeffs + b_coeffs.*d_coeffs)./(c_coeffs.^2 + d_coeffs.^2) ) != 1/2
            @warn "Filter coefficients are not normalised: $(sum((a_coeffs.*c_coeffs + b_coeffs.*d_coeffs)./(c_coeffs.^2 + d_coeffs.^2) )) != 0.5"
        end
    end
    
    # Finally, we can define the grid, if not given (as is typical)
    example_timeseries = FieldTimeSeries(original_data_filename, velocity_names[1]; architecture=architecture, backend=backend)
    grid = isnothing(grid) ? example_timeseries.grid : grid

    filter_mode = "offline"  

    return OfflineFilterConfig(original_data_filename,
                            var_names_to_filter,
                            velocity_names,
                            T_start,
                            T_end,
                            T,
                            architecture,
                            T_out,
                            filter_params,
                            Δt,
                            backend,
                            map_to_mean,
                            forward_output_filename,
                            backward_output_filename,
                            output_filename,
                            npad,
                            delete_intermediate_files,
                            compute_Eulerian_filter,
                            output_netcdf,
                            output_original_data,
                            advection,
                            grid,
                            filter_mode)

end


end # module OfflineLagrangianFilter


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


##### Define a config structure

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
    compute_mean_velocities::Bool
    delete_intermediate_files::Bool
    compute_Eulerian_filter::Bool
    output_netcdf::Bool
    output_original_data::Bool
    advection::Union{AbstractAdvectionScheme, Nothing}
    grid::AbstractGrid
    filter_mode::String
    label::String

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
                        compute_mean_velocities::Bool = true,
                        delete_intermediate_files::Bool = true,
                        compute_Eulerian_filter::Bool = false,
                        output_netcdf::Bool = false,
                        output_original_data::Bool = true,
                        advection::Union{AbstractAdvectionScheme, Nothing} = WENO(),
                        grid::Union{AbstractGrid, Nothing} = nothing,
                        filter_mode::String = "offline",
                        label::String = "")

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
  - `N`, `freq_c`: Parameters for a Butterworth squared filter. `N` is the order of the filter, and `freq_c` is the cutoff frequency. 
     These are used to automatically generate `filter_params` if not provided. Must be specified together if `filter_params` is not given.
  - `filter_params`: A `NamedTuple` containing the coefficients for a custom filter. Only filter_params OR `N` and `freq_c` should be given.
  - `Δt`: The time step for the internal Lagrangian filter simulation. If `nothing`, it defaults to `T_out / 10`, but this may not be appropriate.
  - `backend`: The backend for loading `FieldTimeSeries` data. See `Oceananigans.Fields.FieldTimeSeries`. Default: `InMemory(4)`.
  - `map_to_mean`: A `Bool` indicating whether to map filtered data to the mean position (i.e. calculate generalised Lagrangian mean). Default: `true`.
  - `forward_output_filename`: The filename for the output of the forward filter pass. Default: `"forward_output.jld2"`.
  - `backward_output_filename`: The filename for the output of the backward filter pass. Default: `"backward_output.jld2"`.
  - `output_filename`: The filename for the final combined and mapped output. Default: `"filtered_output.jld2"`.
  - `npad`: The number of cells to pad the interpolation to mean position, used when there are periodic boundary conditions. Default: `5`.
  - `compute_mean_velocities`: A `Bool` indicating whether to compute the mean velocities from the maps. Default: `true`.
  - `delete_intermediate_files`: A `Bool` indicating whether to delete `forward_output.jld2` and `backward_output.jld2` after the final combined file is created. Default: `true`.
  - `compute_Eulerian_filter`: A `Bool` indicating whether to also compute an Eulerian-mean-based filter for comparison. Default: `false`.
  - `output_netcdf`: A `Bool` indicating whether to also convert the final JLD2 output file to a NetCDF file. Default: `false`.
    - `output_original_data`: A `Bool` indicating whether to include the original data in the final output file for comparison. Default: `true`.
  - `advection`: The advection scheme to use for the Lagrangian filter simulation. Default: `WENO()`. Using lower-order schemes may be a source of error.
  - `grid`: The grid for the simulation. If `nothing`, the grid is inferred from the `original_data_filename` (preferred option)
  - `filter_mode`: A `String` indicating whether to run the filter in "offline" or "online" mode. Default: "offline". TODO use multiple dispatch for this instead.
  - `label`: A `String` label for the variables that will be created to pass to the model. For use when multiple filter configurations are to be run
     at the same time.  Default: "".
# Example:

```jldoctest offline config
using OceananigansLagrangianFilter
using Oceananigans.Units
path_to_sim = "../test/data/reference_sim.jld2"
filter_config = OfflineFilterConfig(original_data_filename=path_to_sim, 
                                    output_filename = "output_file.jld2", 
                                    var_names_to_filter = ("b",), 
                                    velocity_names = ("u","w"), 
                                    architecture = CPU(), 
                                    Δt = 20minutes, 
                                    T_out = 1hour, 
                                    N = 2, 
                                    freq_c = 1e-4/2, 
                                    compute_mean_velocities = true, 
                                    output_netcdf = true, 
                                    delete_intermediate_files = true, 
                                    compute_Eulerian_filter = true) 

# output
┌ Info: Advection for Lagrangian filtering will be performed using only velocities ("u", "w") - 
│ any other velocity components will be zero by default. Maps for regridding to mean position will
└ be computed corresponding to velocities: ("u", "w").
[ Info: Mean velocities corresponding to ("u", "w") will be computed.
[ Info: Filter interval will be from T_start=0.0 to T_end=86400.0, duration T=86400.0
[ Info: Setting filter parameters to use Butterworth squared, order 2, cutoff frequency 5.0e-5
OfflineFilterConfig("../test/data/reference_sim.jld2", ("b",), ("u", "w"), 0.0, 86400.0, 86400.0, CPU(), 3600.0, (a1 = 1.767766952966369e-5, b1 = 1.767766952966369e-5, c1 = 3.535533905932738e-5, d1 = 3.535533905932738e-5, N_coeffs = 1), 1200.0, InMemory{Int64}(1, 4), true, "forward_output.jld2", "backward_output.jld2", "output_file.jld2", 5, true, true, true, true, true, WENO{3, Float64, Float32}(order=5)
├── buffer_scheme: WENO{2, Float64, Float32}(order=3)
└── advection_velocity_scheme: Centered(order=4), 10×1×10 RectilinearGrid{Float64, Periodic, Flat, Bounded} on CPU with 3×0×3 halo
├── Periodic x ∈ [-5000.0, 5000.0) regularly spaced with Δx=1000.0
├── Flat y                         
└── Bounded  z ∈ [-100.0, 0.0]     regularly spaced with Δz=10.0, "offline", "")

```

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
                            compute_mean_velocities::Bool = true,
                            delete_intermediate_files::Bool = true,
                            compute_Eulerian_filter::Bool = false,
                            output_netcdf::Bool = false,
                            output_original_data::Bool = true,
                            advection::Union{AbstractAdvectionScheme, Nothing} = WENO(),
                            grid::Union{AbstractGrid, Nothing} = nothing,
                            label::String = "")

    # Check that the original file exists 
    if !isfile(original_data_filename)
        error("Source file not found: $original_data_filename")
    end

    # Check that velocities aren't in the var_names_to_filter
    for vel_str in ("u","v","w")
        if vel_str in var_names_to_filter
            error("Velocity variable '$vel_str' cannot be in 'var_names_to_filter'. 
If you want to filter a velocity that you will also be advecting with, include it in 'velocity_names'. 
If you want to filter a velocity that you don't want to advect with, then output it as a variable with 
a different name in the original simulation.")
        end
    end

    if "t" in var_names_to_filter
        error("Time variable 't' cannot be in 'var_names_to_filter'.")
    end

    # Notify about the velocities that will be used
    if map_to_mean
        @info "Advection for Lagrangian filtering will be performed using only velocities $(velocity_names) - 
any other velocity components will be zero by default. Maps for regridding to mean position will
be computed corresponding to velocities: $(velocity_names)."
    else
        @info "Advection for Lagrangian filtering will be performed using only velocities $(velocity_names) - 
any other velocity components will be zero by default."
    end

    if compute_mean_velocities
        @info "Mean velocities corresponding to $(velocity_names) will be computed."
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
        filter_params = set_offline_BW2_filter_params(;N,freq_c)
        @info "Setting filter parameters to use Butterworth squared, order $N, cutoff frequency $freq_c"
    else # !isnothing(filter_params) && isnothing(N) && isnothing(freq_c)
        # User has specified filter_params directly, but we should check it has the right fields
        if haskey(filter_params, :N_coeffs)
            if filter_params.N_coeffs == 0.5 # Single exponential special case
                if !all((haskey(filter_params, :a1) , haskey(filter_params, :c1)))
                    error("For N_coeffs=0.5, filter_params must have fields :a1, and :c1")
                end
            elseif Int(floor(filter_params.N_coeffs)) !== filter_params.N_coeffs
                error("N_coeffs must be a positive integer or 0.5")
            else

                if !all((haskey(filter_params, Symbol(coeff,i)) for coeff in ["a","b","c","d"] for i in 1:filter_params.N_coeffs))
                    error("For N_coeffs>0.5, filter_params must have fields :a1, :a2, ..., :b1, :b2, ..., :c1, :c2, ..., :d1, :d2, ...")
                end
            
            end
        else # N_coeffs isn't provided, but we might be able to infer it
            if floor(length(filter_params)/4) == length(filter_params)/4
                filter_params = merge(filter_params, (N_coeffs = Int(length(filter_params)/4),))
                # But we still have to check that the right entries are there:
                if !all((haskey(filter_params, Symbol(coeff,i)) for coeff in ["a","b","c","d"] for i in 1:filter_params.N_coeffs))
                    error("filter_params must have fields :a1, :a2, ..., :b1, :b2, ..., :c1, :c2, ..., :d1, :d2, ...")
                end
            elseif length(filter_params) == 2
                filter_params = merge(filter_params, (N_coeffs = 0.5,))
                if !all((haskey(filter_params, :a1) , haskey(filter_params, :c1)))
                    error("For a filter with two coefficients, filter_params must have fields :a1, and :c1")
                end
            else
                error("filter_params must have either 2 entries (for single exponential) or a multiple of 4 entries, e.g. 2*N entries for Butterworth squared of order N, N even.")
            
            end

        end
    end

    # Check normalisation of filter coefficients
    if filter_params.N_coeffs == 0.5
        if !(filter_params.a1*2 ≈ filter_params.c1)
            @warn "Filter coefficients are not normalised: 2*a1=$(2*filter_params.a1) != c1=$(filter_params.c1). 
You can continue, but setting `map_to_mean=false` as the map is now meaningless."
            map_to_mean = false
        end
    else
        a_coeffs = [filter_params[Symbol("a",i)] for i in 1:filter_params.N_coeffs]
        b_coeffs = [filter_params[Symbol("b",i)] for i in 1:filter_params.N_coeffs]
        c_coeffs = [filter_params[Symbol("c",i)] for i in 1:filter_params.N_coeffs] 
        d_coeffs = [filter_params[Symbol("d",i)] for i in 1:filter_params.N_coeffs]
        if !(sum((a_coeffs.*c_coeffs + b_coeffs.*d_coeffs)./(c_coeffs.^2 + d_coeffs.^2) ) ≈ 1/2)
            @warn "Filter coefficients are not normalised: $(sum((a_coeffs.*c_coeffs + b_coeffs.*d_coeffs)./(c_coeffs.^2 + d_coeffs.^2) )) != 0.5
You can continue, but setting `map_to_mean=false` as the map is now meaningless."
            map_to_mean = false
        end
    end

    # Finally, we can define the grid, if not given (as is typical)
    example_timeseries = FieldTimeSeries(original_data_filename, velocity_names[1]; architecture=architecture, backend=backend)
    grid = isnothing(grid) ? example_timeseries.grid : grid

    # Give a warning if the grid has an immersed boundary
    if grid isa ImmersedBoundaryGrid
        @warn "The final interpolation to mean position does not yet work well with immersed boundaries - consider setting map_to_mean=false"
    end

    underlying_rectilinear_grid = (grid isa RectilinearGrid) || ((grid isa ImmersedBoundaryGrid) && (grid.underlying_grid isa RectilinearGrid))

    # Give warning about interpolation if grid is not RectilinearGrid and turn off interpolation for now
    if !underlying_rectilinear_grid && map_to_mean
        @warn "The final interpolation to mean position currently only works for RectilinearGrids - setting map_to_mean=false"
        map_to_mean = false
    end

    underlying_latlon_grid = (grid isa LatitudeLongitudeGrid) || ((grid isa ImmersedBoundaryGrid) && (grid.underlying_grid isa LatitudeLongitudeGrid))

    if !underlying_rectilinear_grid && !underlying_latlon_grid && output_netcdf
        @warn "The grid type $(typeof(grid)) won't work with NetCDF output functionality - setting output_netcdf=false"
        output_netcdf = false
    end

    # Make sure that map_to_mean is false if advection is nothing (Eulerian filter)

    if isnothing(advection) && map_to_mean
        @warn "Advection scheme is 'nothing' (Eulerian filter) so setting map_to_mean=false"
        map_to_mean = false
    elseif isnothing(advection) 
        @info "Advection scheme is 'nothing' so the Eulerian (not Lagrangian) filter will be computed."
    end

    # Warn if Eulerian filter is being calculated twice
    if compute_Eulerian_filter && isnothing(advection)
        @warn "compute_Eulerian_filter=true and advection is 'nothing' - Eulerian filter will be computed twice, so we'll set compute_Eulerian_filter=false."
        compute_Eulerian_filter = false
    end

    
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
                            compute_mean_velocities,
                            delete_intermediate_files,
                            compute_Eulerian_filter,
                            output_netcdf,
                            output_original_data,
                            advection,
                            grid,
                            filter_mode,
                            label)

end


end # module OfflineLagrangianFilter


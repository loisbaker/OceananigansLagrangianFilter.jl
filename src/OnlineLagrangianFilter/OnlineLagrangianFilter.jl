module OnlineLagrangianFilter

using ..OceananigansLagrangianFilter: AbstractConfig
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Architectures

using ..Utils

export OnlineFilterConfig
"""
    OnlineFilterConfig(;

    )

A configuration object for online filtering.
"""
struct OnlineFilterConfig <: AbstractConfig
    grid::AbstractGrid
    output_filename::String
    var_names_to_filter::Tuple{Vararg{String}}
    velocity_names::Tuple{Vararg{String}} 
    filter_params::NamedTuple
    map_to_mean::Bool
    compute_mean_velocities::Bool
    npad::Int
    filter_mode::String
    
end

"""
    OnlineFilterConfig(; grid::AbstractGrid,
                            output_filename::String = "online_filtered_output",
                            var_names_to_filter::Tuple{Vararg{String}},
                            velocity_names::Tuple{Vararg{String}},
                            N::Union{Int, Nothing} = nothing,
                            freq_c::Union{Int, Nothing} = nothing,
                            filter_params::Union{NamedTuple, Nothing} = nothing,
                            map_to_mean::Bool = true,
                            compute_mean_velocities::Bool = true,
                            npad::Int = 5,
                            filter_mode::String = "online"
                            )

Constructs a configuration object for offline Lagrangian filtering of Oceananigans data.
This function validates the input data file, time specifications, and filter parameters
before creating the `OfflineFilterConfig` object.

Keyword arguments
=================
  - `grid`: (required) The grid for the simulation. If `nothing`, the grid is inferred from the `original_data_filename` (preferred option)
  - `output_filename`: The filename for the output of the online filtered data. Default: `"online_filtered_output"`.
  - `var_names_to_filter`: (required) A `Tuple` of `String`s specifying the names of the tracer variables to be filtered.
  - `velocity_names`: (required) A `Tuple` of `String`s specifying the names of the velocity fields in the data file to be used for advection.
  - `N`, `freq_c`: Parameters for a Butterworth filter. `N` is the order of the filter, and `freq_c` is the cutoff frequency. 
     These are used to automatically generate `filter_params` if not provided. Must be specified together if `filter_params` is not given.
  - `filter_params`: A `NamedTuple` containing the coefficients for a custom filter. Only filter_params OR `N` and `freq_c` should be given.
  - `map_to_mean`: A `Bool` indicating whether to map filtered data to the mean position (i.e. calculate generalised Lagrangian mean). Default: `true`.
  - `compute_mean_velocities`: A `Bool` indicating whether to compute the mean velocities from the maps. Default: `true`.
  - `npad`: The number of cells to pad the interpolation to mean position, used when there are periodic boundary conditions. Default: `5`.
  - `compute_Eulerian_filter`: A `Bool` indicating whether to also compute an Eulerian-mean-based filter for comparison. Default: `false`.
- `filter_mode`: A `String` indicating whether to run the filter in "offline" or "online" mode. Default: "online". TODO use multiple dispatch for this instead.
"""
function OnlineFilterConfig(; grid::AbstractGrid,
                            output_filename::String = "online_filtered_output",
                            var_names_to_filter::Tuple{Vararg{String}},
                            velocity_names::Tuple{Vararg{String}},
                            N::Union{Int, Nothing} = nothing,
                            freq_c::Union{Real, Nothing} = nothing,
                            filter_params::Union{NamedTuple, Nothing} = nothing,
                            map_to_mean::Bool = true,
                            compute_mean_velocities::Bool = true,
                            npad::Int = 5,
                            )

    # Make sure we have some filter parameters
    if !isnothing(filter_params) && (!isnothing(N) || !isnothing(freq_c))
        error("Specify either filter_params or N and freq_c, not both.")
    elseif isnothing(filter_params) && (isnothing(N) || isnothing(freq_c))
        error("Must specify either filter_params or both N and freq_c.")
    elseif isnothing(filter_params) && !isnothing(N) && !isnothing(freq_c)
        filter_params = set_online_BW_filter_params(;N,freq_c)
        @info "Setting filter parameters to use Butterworth order $N, cutoff frequency $freq_c"
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
                error("filter_params must have either 2 entries (for single exponential) or a multiple of 4 entries, e.g. 2*N entries for Butterworth squared of order N, N even.")
            
            end

        end
    end

    # Check normalisation of filter coefficients
    if filter_params.N_coeffs == 0.5
        if !(filter_params.a1 ≈ filter_params.c1)
            @warn "Filter coefficients are not normalised: a1=$(filter_params.a1) != c1=$(filter_params.c1)"
        end
    else
        a_coeffs = [filter_params[Symbol("a",i)] for i in 1:filter_params.N_coeffs]
        b_coeffs = [filter_params[Symbol("b",i)] for i in 1:filter_params.N_coeffs]
        c_coeffs = [filter_params[Symbol("c",i)] for i in 1:filter_params.N_coeffs] 
        d_coeffs = [filter_params[Symbol("d",i)] for i in 1:filter_params.N_coeffs]
        if !(sum((a_coeffs.*c_coeffs + b_coeffs.*d_coeffs)./(c_coeffs.^2 + d_coeffs.^2) ) ≈ 1.0)
            @warn "Filter coefficients are not normalised: $(sum((a_coeffs.*c_coeffs + b_coeffs.*d_coeffs)./(c_coeffs.^2 + d_coeffs.^2) )) != 1"
        end
    end

    filter_mode = "online"
    return OnlineFilterConfig(grid,
                            output_filename,
                            var_names_to_filter,
                            velocity_names,
                            filter_params,
                            map_to_mean,
                            compute_mean_velocities,
                            npad,
                            filter_mode
                            )
 
end


end # module OnlineLagrangianFilter



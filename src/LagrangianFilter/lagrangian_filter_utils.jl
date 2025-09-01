"""
    copy_file_metadata!(original_file::JLD2.JLDFile, new_file::JLD2.JLDFile, timeseries_vars_to_copy::Tuple{Vararg{String}})

Copies essential metadata from an existing Oceananigans JLD2 output file to a new file.

This function is a utility for creating a new file structure with all necessary
simulation metadata (such as grid information and serialized objects) without
copying large timeseries data. It copies core simulation metadata (`grid`, `serialized`,
`coriolis`) and the `serialized` entries for specified timeseries variables.

# Arguments
- `original_file::JLD2.JLDFile`: The source JLD2 file.
- `new_file::JLD2.JLDFile`: The destination JLD2 file.
- `timeseries_vars_to_copy::Tuple{Vararg{String}}`: A tuple of timeseries variable names e.g., ("u", "v").

"""
function copy_file_metadata!(original_file::JLD2.JLDFile, new_file::JLD2.JLDFile, timeseries_vars_to_copy::Tuple{Vararg{String}})

    # Copy over metadata that isn't associated with variables
    for path in ("grid","serialized","coriolis")
        _copy_jld2_recursive!(original_file, new_file, path)
    end

    # Copy serialized entries for timeseries variables
    for var in timeseries_vars_to_copy
        _copy_jld2_recursive!(original_file, new_file, "timeseries/$var/serialized")
    end

end



"""
    _copy_jld2_recursive!(source::JLD2.JLDFile, dest::JLD2.JLDFile, path::String)

A helper function to recursively copy data and groups.
"""
function _copy_jld2_recursive!(source::JLD2.JLDFile, dest::JLD2.JLDFile, path::String)
    
    # Check if the path points to a JLD2 group
    if isa(source[path], JLD2.Group)
        # Recreate the group in the destination file
        dest_group = JLD2.Group(dest, path)

        # Recursively copy the contents of the group
        for key in keys(source[path])
            _copy_jld2_recursive!(source, dest, "$path/$key")
        end
    else
        # Simply copy the variable if it's not a group
        dest[path] = source[path]
    end
end

"""
    create_input_data_on_disk(original_data_filename::String,
                              var_names_to_filter::Tuple{Vararg{String}},
                              velocity_names::Tuple{Vararg{String}};
                              direction::String="forward",
                              T_start::Union{Real, Nothing}=nothing,
                              T_end::Union{Real, Nothing}=nothing)

Creates a new JLD2 file for Lagrangian filtering from an existing Oceananigans
simulation output file by subsetting and re-ordering data on disk.

This function copies metadata and a truncated time-series of specified variables
into a new file, `<original_data_filename>_filter_input.jld2`. The time coordinate
and velocity fields are transformed based on `direction`.

# Functionality
- **Data Subsetting:** Copies specified variables over a time range `[T_start, T_end]`.
- **Direction Re-ordering:**
    - `"forward"`: Time is shifted to start at `t=0`.
    - `"backward"`: Time is inverted (`t = T_end - t_simulation`) and data is reversed,
      with velocities negated.

# Arguments
- `original_data_filename::String`: Path to the source JLD2 file.
- `var_names_to_filter::Tuple{Vararg{String}}`: Names of tracer variables to copy e.g., ("T", "S").
- `velocity_names::Tuple{Vararg{String}}`: Names of velocity variables to copy e.g., ("u", "v", "w").

# Keyword Arguments
- `direction::String`: The direction of the new time coordinate (`"forward"` or `"backward"`).
- `T_start::Union{Real, Nothing}`: Start time for the data subset. Defaults to the first time step.
- `T_end::Union{Real, Nothing}`: End time for the data subset. Defaults to the last time step.

# Returns
- `T_filter::Real`: The total time span of the filtered data (`T_end - T_start`).

# Example
```julia
# Create a new JLD2 file with filtered `u` and `v` velocities and a `T` tracer
# over a specific time period, with the time coordinate running backwards.
T_filter = create_input_data_on_disk(
    "simulation.jld2",
    ("T",),
    ("u", "v");
    direction="backward",
    T_start=100.0,
    T_end=200.0
)
```
"""
function create_input_data_on_disk(original_data_filename::String, 
                                   var_names_to_filter::Tuple{Vararg{String}}, 
                                   velocity_names::Tuple{Vararg{String}}; 
                                   direction::String="forward", 
                                   T_start::Union{Real, Nothing}=nothing, 
                                   T_end::Union{Real, Nothing}=nothing)
    
    # Check that a valid direction has been given
    if direction ∉ ("backward", "forward")
        error("Invalid direction: $direction. Must be 'backward' or 'forward'.")
    end

    # Check that the original file exists 
    if !isfile(original_data_filename)
        error("Source file not found: $original_data_filename")
    end

    # Create a new filename and delete any file that already has that name
    new_filename = original_data_filename[1:end-5]*"_filter_input.jld2"
    if isfile(new_filename)
        rm(new_filename)
    end

    # Open the original file for reading 
    jldopen(original_data_filename,"r") do original_file

        # Check if velocities and tracers given are in the original_file
        for var in (var_names_to_filter..., velocity_names...)
            if !haskey(original_file["timeseries"], var)
                error("Variable '$var' not found in original data file.")
            end
        end

        # Check if T_start and T_end are found in the original_file
        iterations = parse.(Int, keys(original_file["timeseries/t"]))
        times = [original_file["timeseries/t/$iter"] for iter in iterations]
        
        # If T_start is not specified, use the first simulation time
        T_start = isnothing(T_start) ? times[1] : T_start

        # If T_end is not specified, use the last simulation time
        T_end = isnothing(T_end) ? times[end] : T_end

        if T_start < times[1] || T_start > times[end]
            error("T_start=$T_start is outside the range of the original data: [$times[1], $times[end]].")
        end

        if T_end < times[1] || T_end > times[end]
            error("T_end=$T_end is outside the range of the original data: [$times[1], $times[end]].")
        end

        # Create a truncated iterations variable with just the times to copy
        iterations_truncated = iterations[(times .>= T_start) .& (times .<= T_end)]

        # Create a new filename for the filtered input data
        jldopen(new_filename, "w") do new_file

            # We need to copy the standard metadata from the old file to the new file
            copy_file_metadata!(original_file, new_file, (var_names_to_filter..., velocity_names...))

            # Add an indicator of the direction 
            new_file["direction"] = direction

            # Add an old (t_simulation) and a new (t) time variable 
            t_sim_group = JLD2.Group(new_file, "timeseries/t_simulation")
            t_group = JLD2.Group(new_file, "timeseries/t")

            if direction == "forward"

                # First write times, starting from T_start
                for iter in iterations_truncated
                    t_sim_group["$iter"] = original_file["timeseries/t/$iter"]
                    t_group["$iter"] = original_file["timeseries/t/$iter"] - T_start
                end

                # Then write in data, as in original file
                for var in (var_names_to_filter..., velocity_names...)
                    for iter in iterations_truncated
                        new_file["timeseries/$var/$iter"] = original_file["timeseries/$var/$iter"]
                    end
                end

            elseif direction == "backward"

                # First write times, starting from T_end
                for iter in reverse(iterations_truncated)
                    t_sim_group["$iter"] = original_file["timeseries/t/$iter"]
                    t_group["$iter"] = T_end - original_file["timeseries/t/$iter"]
                end

                # Then write in data (as in original file for tracers and negated for velocities)
                for var in var_names_to_filter
                    for iter in reverse(iterations_truncated)
                        new_file["timeseries/$var/$iter"] = original_file["timeseries/$var/$iter"]
                    end
                end

                for var in velocity_names
                    for iter in reverse(iterations_truncated)
                        new_file["timeseries/$var/$iter"] = -original_file["timeseries/$var/$iter"]
                    end
                end

                
            end

            # Calculate and return the total filter time
            T_filter = T_end - T_start
            return T_filter
        end
    end
end

"""
    load_data(original_data_filename::String,
              var_names_to_filter::Tuple{Vararg{String}},
              velocity_names::Tuple{Vararg{String}};
              architecture=CPU(),
              backend=InMemory())

Loads time-series data for the Lagrangian filter from a JLD2 file.

This function assumes that the input data has already been prepared by
`create_input_data_on_disk`. It loads the specified `FieldTimeSeries`
and extracts the grid from the first loaded time series.

# Arguments
- `original_data_filename::String`: The path to the original JLD2 file.
- `var_names_to_filter::Tuple{Vararg{String}}`: A tuple of tracer variable names.
- `velocity_names::Tuple{Vararg{String}}`: A tuple of velocity variable names.

# Keyword Arguments
- `architecture`: The architecture on which to load the data (`CPU` or `GPU`).
  Defaults to `CPU`.
- `backend`: The data backend for the `FieldTimeSeries` (`InMemory` or `OnDisk`).
  Defaults to `InMemory`.

# Returns
- `Tuple{Tuple, Tuple, Any}`: A tuple containing:
    - `velocity_timeseries::Tuple`: A tuple of loaded velocity `FieldTimeSeries`.
    - `var_timeseries::Tuple`: A tuple of loaded tracer `FieldTimeSeries`.
    - `grid`: The grid object from the loaded data.
"""
function load_data(original_data_filename::String, 
                   var_names_to_filter::Tuple{Vararg{String}}, 
                   velocity_names::Tuple{Vararg{String}}; 
                   architecture=CPU(), 
                   backend=InMemory())

    input_data_filename = original_data_filename[1:end-5] * "_filter_input.jld2"
    
    velocity_timeseries = Tuple(FieldTimeSeries(input_data_filename, name; architecture=architecture, backend=backend) for name in velocity_names)
    var_timeseries = Tuple(FieldTimeSeries(input_data_filename, name; architecture=architecture, backend=backend) for name in var_names_to_filter)

    grid = velocity_timeseries[1].grid

    return velocity_timeseries, var_timeseries, grid
end

"""
    set_BW_filter_params(;N::Int=1, freq_c::Real=1)

Calculates the coefficients for a filter with that has frequency response given by a
Butterworth filter with order `2^N` and cutoff frequency `freq_c`, squared. 

Uses 2^N exponentials and 2^(N-1) sets of coefficients (which come in pairs to ensure real equations)

Frequency response: Ghat(omega) = 1 / (1 + (omega / freq_c)^(2^(N+1)))
Real filter shape: G(t) = sum_{i=1}^{2^(N-1)} exp(-c_i*abs(t))*(a_i*cos(d_i * abs(t)) + b_i*sin(d_i * abs(t)))

# Keyword Arguments
- `N::Int`: log_2 order of the Butterworth filter. Must be a positive integer.
  Defaults to 1.
- `freq_c::Real`: The cutoff frequency of the filter. Defaults to 1.

# Returns
- `NamedTuple`: A `NamedTuple` containing the calculated coefficients.
  The coefficients are named `a1, b1, c1, d1, a2, b2, ...` up to `2^(N-1)` pairs.

"""
function set_BW_filter_params(;N::Int=1,freq_c::Real=1) 
    N_coeffs = 2^(N-1)
    filter_params = NamedTuple()
    for i in 1:N_coeffs
        
        a = (freq_c/2^N)*sin(pi/(2^(N+1))*(2*i-1))
        b = (freq_c/2^N)*cos(pi/(2^(N+1))*(2*i-1))
        c = freq_c*sin(pi/(2^(N+1))*(2*i-1))
        d = freq_c*cos(pi/(2^(N+1))*(2*i-1))

        temp_params = NamedTuple{(Symbol("a$i"), Symbol("b$i"),Symbol("c$i"),Symbol("d$i"))}([a,b,c,d])
        filter_params = merge(filter_params,temp_params)
    end
    
    return filter_params
end

"""
    create_original_vars(var_names_to_filter::Tuple{Vararg{String}}, grid::AbstractGrid)

Creates a `NamedTuple` of `CenterField`s on the given `grid` for each variable name.

# Arguments
- `var_names_to_filter::Tuple{Vararg{String}}`: A tuple of variable names e.g., ("T", "S").
- `grid`: The original simulation grid 

# Returns
- `NamedTuple`: A `NamedTuple` with keys corresponding to `var_names_to_filter`
  and values as `CenterField` objects initialized on `grid`.

"""
function create_original_vars(var_names_to_filter::Tuple{Vararg{String}}, grid::AbstractGrid)
    # Creates auxiliary fields to store the saved variables
    vars = Dict()
    for var_name in var_names_to_filter
        vars[Symbol(var_name)] = CenterField(grid)
    end
    return NamedTuple(vars)
end

"""
    create_filtered_vars(var_names_to_filter::Tuple{Vararg{String}},
                         velocity_names::Tuple{Vararg{String}},
                         filter_params::NamedTuple;
                         map_to_mean::Bool=true)

Creates a tuple of `Symbol`s for the filtered variables used in the
Lagrangian filter. The symbols are generated to represent the
cosine (`_C#`) and sine (`_S#`) components of the final filtered variables.

The number of components for each variable is determined by the number of
coefficients in `filter_params`. If `map_to_mean` is true, additional symbols
are created for map variables to interpolate to mean position.

# Arguments
- `var_names_to_filter::Tuple{Vararg{String}}`: Names of tracer variables to filter.
- `velocity_names::Tuple{Vararg{String}}`: Names of velocity variables.
- `filter_params::NamedTuple`: A `NamedTuple` of filter coefficients, typically
  from `set_BW_filter_params`, or custom.

# Keyword Arguments
- `map_to_mean::Bool`: If `true`, includes symbols for map variables.
  Defaults to `true`.

# Returns
- `Tuple{Vararg{Symbol}}`: A tuple of symbols representing the filtered variables.
  The symbols for the cosine components are listed first, followed by the sine components.

"""
function create_filtered_vars(var_names_to_filter::Tuple{Vararg{String}},
                              velocity_names::Tuple{Vararg{String}},
                              filter_params::NamedTuple;
                              map_to_mean::Bool=true)
                              
    N_coeffs = Int(length(filter_params) / 4)

    gC_symbols = Symbol[]
    gS_symbols = Symbol[]

    for var_name in var_names_to_filter
        for i in 1:N_coeffs
            push!(gC_symbols, Symbol(var_name, "C", i))
            push!(gS_symbols, Symbol(var_name, "S", i))
        end
    end

    # May also need xi maps. We need one for every velocity dimension, so lets use the velocity names to name them
    if map_to_mean
        for vel_name in velocity_names
            for i in 1:N_coeffs
                push!(gC_symbols, Symbol("xi_", vel_name, "_C", i))
                push!(gS_symbols, Symbol("xi_", vel_name, "_S", i))
            end
        end
    end

    return (Tuple(gC_symbols)..., Tuple(gS_symbols)...)
end

"""
    make_gC_forcing_func(i::Int)

Creates a forcing function for the `i`-th cosine component of a filtered variable.
The returned function can be passed directly to `ModelForcing` in Oceananigans.

# Arguments
- `i::Int`: The index of the filter coefficient pair to use.

# Returns
- A function with the signature `(x, y, z, t, gC, gS, p)` (in 3D) that returns the forcing
  function for `gCi`.
"""
function make_gC_forcing_func(i::Int)
    
    function gC_forcing(args...)
        filter_params = args[end]
        c = getproperty(filter_params,Symbol("c$i"))
        d = getproperty(filter_params,Symbol("d$i"))
        gC = args[end-2]
        gS = args[end-1]
        return -c .* gC .- d .* gS
    end
    return gC_forcing
end

"""
    make_gS_forcing_func(i::Int)

Creates a forcing function for the `i`-th sine component of a filtered variable.
The returned function can be passed directly to `ModelForcing` in Oceananigans.

# Arguments
- `i::Int`: The index of the filter coefficient pair to use.

# Returns
- A function with the signature `(x, y, z, t, gC, gS, p)` (in 3D) that returns the forcing
  function for `gSi`.
"""
function make_gS_forcing_func(i::Int)
    function gS_forcing(args...)
        filter_params = args[end]
        c = getproperty(filter_params,Symbol("c$i"))
        d = getproperty(filter_params,Symbol("d$i"))
        gC = args[end-2]
        gS = args[end-1]
        return -c .* gS .+ d .* gC
    end
    return gS_forcing
end

"""
    make_xiC_forcing_func(i::Int)

Creates a forcing function for the `i`-th cosine component of a mapping variable.
The returned function can be passed directly to `ModelForcing` in Oceananigans.

# Arguments
- `i::Int`: The index of the filter coefficient pair to use.

# Returns
- A function with the signature `(x, y, z, t, u, p)` (in 3D) that returns the forcing
  function for `xi_u_Ci`, where u is the corresponding velocity.
"""
function make_xiC_forcing_func(i::Int)
    function xiC_forcing(args...)
        filter_params = args[end]
        c = getproperty(filter_params,Symbol("c$i"))
        d = getproperty(filter_params,Symbol("d$i"))
        u = args[end-1]

        return -c .* u ./ (c.^2 .+ d.^2)
    end
    return xiC_forcing
end

"""
    make_xiS_forcing_func(i::Int)

Creates a forcing function for the `i`-th sine component of a mapping variable.
The returned function can be passed directly to `ModelForcing` in Oceananigans.

# Arguments
- `i::Int`: The index of the filter coefficient pair to use.

# Returns
- A function with the signature `(x, y, z, t, u, p)` (in 3D) that returns the forcing
  function for `xi_u_Si`, where u is the corresponding velocity
"""
function make_xiS_forcing_func(i::Int)
    function xiS_forcing(args...)
        filter_params = args[end]
        c = getproperty(filter_params,Symbol("c$i"))
        d = getproperty(filter_params,Symbol("d$i"))
        u = args[end-1]
        
        return -d .* u ./ (c.^2 .+ d.^2)
    end
    return xiS_forcing
end

"""
    create_forcing(filtered_vars::Tuple{Vararg{Symbol}},
                   var_names_to_filter::Tuple{Vararg{String}},
                   velocity_names::Tuple{Vararg{String}},
                   filter_params::NamedTuple)

Creates a `NamedTuple` of `Forcing` objects for each filtered variable
to be used in an Oceananigans model.

# Arguments
- `filtered_vars::Tuple{Vararg{Symbol}}`: A tuple of symbols for all
  filtered variables, from `create_filtered_vars`.
- `var_names_to_filter::Tuple{Vararg{String}}`: Names of tracer variables
  to be filtered (user defined).
- `velocity_names::Tuple{Vararg{String}}`: Names of velocity variables
  to be filtered (user defined).
- `filter_params::NamedTuple`: A `NamedTuple` of filter coefficients.

# Returns
- `forcing::NamedTuple`: A `NamedTuple` where each key is a filtered
  variable symbol and each value is an Oceananigans `Forcing` object.
  The forcings are multi-term forcings where appropriate.
"""
function create_forcing(filtered_vars::Tuple{Vararg{Symbol}},
                        var_names_to_filter::Tuple{Vararg{String}},
                        velocity_names::Tuple{Vararg{String}},
                        filter_params::NamedTuple)

    N_coeffs = Int(length(filter_params)/4)

    # Initialize dictionary
    gC_forcings_dict = Dict()
    gS_forcings_dict = Dict()

    # Make a simple forcing function for original data forcing - the final argument is the field dependence.
    original_var_forcing_func(args...) = args[end]

    for var_name in var_names_to_filter
        original_var = Symbol(var_name)
        for i in 1:N_coeffs
            
            gCkey = Symbol(var_name,"C",i)   
            gSkey = Symbol(var_name,"S",i)   
            
            # The forcing for gC is the sum of a filter forcing term and the original data forcing
            gC_forcing_i = Forcing(make_gC_forcing_func(i), parameters = filter_params, field_dependencies = (gCkey,gSkey))
            gC_original_var_forcing = Forcing(original_var_forcing_func,field_dependencies= (;original_var))
            gC_forcings_dict[gCkey] = (gC_forcing_i, gC_original_var_forcing)

            # The forcing for gS is just a filter forcing term
            gS_forcing_i = Forcing(make_gS_forcing_func(i), parameters = filter_params, field_dependencies = (gCkey,gSkey))
            gS_forcings_dict[gSkey] = gS_forcing_i

        end
    end

    # Check if we need xi forcing (implied by number of filtered_vars)
    if length(filtered_vars) > N_coeffs*length(var_names_to_filter)*2
        # Create forcing for xi maps (one for each velocity component)
        for vel_name in velocity_names
            vel_key = Symbol(velocity_name)
            for i in 1:N_coeffs
                gCkey = Symbol("xi_", vel_name, "_C", i) 
                gSkey = Symbol("xi_", vel_name, "_S", i) 
                
                # The forcing for xiC includes a term involving xiC and xiS (as for tracers, so reuse forcing constructor) and a 
                # term involving the corresponding velocity
                gC_forcing_i = Forcing(make_gC_forcing_func(i), parameters =filter_params, field_dependencies = (gCkey,gSkey))
                xiC_forcing_i = Forcing(make_xiC_forcing_func(i),parameters =filter_params, field_dependencies = (;vel_key))
                gC_forcings_dict[gCkey] = (xiC_forcing_i, gC_forcing_i)

                # The forcing for xiS also includes a term involving xiC and xiS and a term involving the corresponding velocity
                gS_forcing_i = Forcing(make_gS_forcing_func(i), parameters =filter_params, field_dependencies = (gCkey,gSkey))
                xiS_forcing_i = Forcing(make_xiS_forcing_func(i),parameters =filter_params, field_dependencies = (;vel_key))
                gSdict[gSkey] = (xiS_forcing_i,gS_forcing_i)

            end
        end
    end
    forcing = (; NamedTuple(gC_forcings_dict)..., NamedTuple(gS_forcings_dict)...)

    return forcing
end

"""
    create_output_fields(model::LagrangianFilter,
                         var_names_to_filter::Tuple{Vararg{String}},
                         velocity_names::Tuple{Vararg{String}},
                         filter_params::NamedTuple)

Creates a `Dict` of the final, physically meaningful output fields by combining
the filtered components (`_C#`, `_S#`) according to the filter coefficients.

# Arguments
- `model`: The Oceananigans Lagrangian Filter model containing the filtered fields.
- `var_names_to_filter::Tuple{Vararg{String}}`: The names of the tracer
  variables that were filtered.
- `velocity_names::Tuple{Vararg{String}}`: The names of the velocity
  variables that were filtered.
- `filter_params::NamedTuple`: The filter coefficients used in the
  model run, typically from `set_BW_filter_params`.

# Returns
- `Dict`: A dictionary mapping descriptive names (e.g., `"T_filtered"`, `"u_filtered"`)
  to the reconstructed output fields. This dictionary also includes the
  original fields for comparison.
"""
function create_output_fields(model::LagrangianFilter,
                              var_names_to_filter::Tuple{Vararg{String}},
                              velocity_names::Tuple{Vararg{String}},
                              filter_params::NamedTuple)

    N_coeffs = Int(length(filter_params) / 4)
    N_filtered_vars = length(model.tracers)
    outputs_dict = Dict()

    for original_var in var_names_to_filter
        
        # Reconstruct the filtered tracer fields, starting with the first coefficient
        gC1 = getproperty(model.tracers,Symbol(original_var * "C1"))
        gS1 = getproperty(model.tracers,Symbol(original_var * "S1"))
        g_total = filter_params.a1 .* gC1 .+ filter_params.b1 .* gS1
    
        # Then add the other coefficients
        for i in 2:N_coeffs
            a = getproperty(filter_params,Symbol("a$i"))
            b = getproperty(filter_params,Symbol("b$i"))
            gCi = getproperty(model.tracers,Symbol(original_var * "C" * i))
            gSi = getproperty(model.tracers,Symbol(original_var * "S" * i))
            g_total .+= a .* gCi .+ b .* gSi
        end
        outputs_dict[original_var * "_filtered"] = g_total
    end

    # Reconstruct the maps, if they exist
    if N_filtered_vars > N_coeffs*length(var_names_to_filter)*2
        for vel in velocity_names

            # Start with the first coefficient
            xiC1 = getproperty(model.tracers,Symbol("xi_" * vel * "_C1"))
            xiS1 = getproperty(model.tracers,Symbol("xi_" * vel * "_S1"))
            g_total = filter_params.a1 .* xiC1 .+ filter_params.b1 .* xiS1

            # Then add the other coefficients
            for i in 2:N_coeffs
                a = getproperty(filter_params,Symbol("a$i"))
                b = getproperty(filter_params,Symbol("b$i"))
                xiCi = getproperty(model.tracers,Symbol("xi_" * vel * "C" * i))
                xiSi = getproperty(model.tracers,Symbol("xi_" * vel * "S" * i))
                g_total .+= a .* xiCi .+ b .* xiSi
            end
            outputs_dict["xi_" * vel] = g_total
        end
    end

    # Let's also add the saved vars for comparison
    for original_var in var_names_to_filter
        outputs_dict[original_var] = getproperty(model.auxiliary_fields, Symbol(original_var))
    end

    return outputs_dict
end

"""
    update_input_data!(sim::Simulation, input_data::NamedTuple)

Updates the input data fields in the Oceananigans model at the current simulation
time. This function is used as a `callback` to ensure the model's velocity and
auxiliary fields are kept in sync with the time-dependent input data.

# Arguments
- `sim::Simulation`: The Oceananigans `Simulation` object.
- `input_data::NamedTuple`: A `NamedTuple` containing the time series data for
  velocities and original tracer variables.

# Side Effects
- Modifies the fields within `sim.model.velocities` and `sim.model.auxiliary_fields`
  in place, setting them to the values from `input_data` at the current simulation time.
"""
function update_input_data!(sim::Simulation, input_data::NamedTuple)
    velocity_timeseries = input_data.velocities
    original_var_timeseries = input_data.original_vars
    model = sim.model
    t = sim.model.clock.time
    
    # Update the velocities
    kwargs = (; (Symbol(vel_fts.name) => vel_fts[Time(t)] for vel_fts in velocity_timeseries)...)
    set!(model; kwargs...)          

    # We also update the saved original variables to be used for forcing - these are auxiliary fields so need to be set separately
    for original_var_fts in original_var_timeseries
        set!(getproperty(model.auxiliary_fields, Symbol(original_var_fts.name)), original_var_fts[Time(t)])
        # halo regions get filled automatically
    end
end

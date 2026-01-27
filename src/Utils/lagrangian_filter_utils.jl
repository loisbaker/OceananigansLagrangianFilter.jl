using Oceananigans: AbstractModel

"""
    copy_file_metadata!(original_file::JLD2.JLDFile, new_file::JLD2.JLDFile,
                        timeseries_vars_to_copy::Tuple{Vararg{String}})

Copies essential metadata from an existing Oceananigans JLD2 output file to a
new file.

This function is a utility for creating a new file structure with all
necessary simulation metadata (such as grid information and serialized
objects) without copying large timeseries data. It copies core simulation
metadata (`grid`, `serialized`, `coriolis`) and the `serialized` entries for
specified timeseries variables.

# Arguments
- `original_file::JLD2.JLDFile`: The source JLD2 file.
- `new_file::JLD2.JLDFile`: The destination JLD2 file.
- `timeseries_vars_to_copy::Tuple{Vararg{String}}`: A tuple of timeseries
  variable names e.g., ("u", "v").

"""
function copy_file_metadata!(original_file::JLD2.JLDFile, new_file::JLD2.JLDFile, 
                             timeseries_vars_to_copy::Tuple{Vararg{String}})

    # Copy over metadata that isn't associated with variables
    for path in ("grid","serialized")
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
    create_input_data_on_disk(config::AbstractConfig; direction::String="forward")

Prepares a new JLD2 file on disk with the time-filtered data for a forward or
backward Lagrangian simulation.

This function performs the following steps:
1.  **Validates `direction`**: Ensures the direction is either `"forward"`
    or `"backward"`.
2.  **Creates a new file**: A new JLD2 file is created with the suffix
    `_filter_input.jld2` and any existing file with the same name is deleted.
3.  **Copies metadata**: Key metadata from the original file (e.g., `grid`
    information) is copied to the new file to maintain consistency.
4.  **Time truncation**: The data is truncated to the time range specified by
    `config.T_start` and `config.T_end`.
5.  **Time shifting**:
    -   For `"forward"` filtering, a new time variable is created, shifted so
        that `t=0` corresponds to `config.T_start`.
    -   For `"backward"` filtering, the data is re-ordered and a new time
        variable is created, shifted so that `t=0` corresponds to `config.T_end`.
6.  **Velocity reversal**: For `"backward"` filtering, the velocity fields
    (`u`, `v`, `w`) are negated to correctly simulate backward advection.

Arguments
=========

- `config`: An instance of `AbstractConfig` containing the file paths,
  variable names, and time specifications.

Keyword Arguments
=================

- `direction`: A `String` indicating the simulation direction. It must be
  either `"forward"` (the default) or `"backward"`.
"""
function create_input_data_on_disk(config::AbstractConfig; direction::String="forward")
    
    original_data_filename = config.original_data_filename
    var_names_to_filter = config.var_names_to_filter
    velocity_names = config.velocity_names
    T_start = config.T_start
    T_end = config.T_end

    # Check that a valid direction has been given
    if direction ∉ ("backward", "forward")
        error("Invalid direction: $direction. Must be 'backward' or 'forward'.")
    end

    # Create a new filename and delete any file that already has that name
    new_filename = original_data_filename[1:end-5]*"_filter_input.jld2"
    if isfile(new_filename)
        rm(new_filename)
    end

    # Open the original file for reading 
    jldopen(original_data_filename,"r") do original_file

        # Check if T_start and T_end are found in the original_file
        iterations = parse.(Int, keys(original_file["timeseries/t"]))
        times = [original_file["timeseries/t/$iter"] for iter in iterations]
        
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

        end
    end
end

"""
    load_data(config::AbstractConfig)

Loads the velocity and tracer data from the intermediate input file created by
`create_input_data_on_disk`. The data for each variable is loaded as a `FieldTimeSeries`
and returned as a single `NamedTuple`.

Arguments
=========

- `config`: An instance of `AbstractConfig` containing the file path, variable names, 
architecture, and backend.

Returns
=======

A `NamedTuple` with fields `velocity_data` and `var_data`, where each field
contains a `Tuple` of `FieldTimeSeries` objects.
"""
function load_data(config::AbstractConfig)

    original_data_filename = config.original_data_filename
    var_names_to_filter = config.var_names_to_filter
    velocity_names = config.velocity_names
    architecture = config.architecture
    backend = config.backend

    input_data_filename = original_data_filename[1:end-5] * "_filter_input.jld2"
    
    velocity_timeseries = Tuple(FieldTimeSeries(input_data_filename, name; architecture=architecture, backend=backend) for name in velocity_names)
    var_timeseries = Tuple(FieldTimeSeries(input_data_filename, name; architecture=architecture, backend=backend) for name in var_names_to_filter)
    input_data = (velocity_data = velocity_timeseries, var_data = var_timeseries)
    return input_data
end


"""
    set_offline_BW2_filter_params(; N::Int=1, freq_c::Real=1)

Calculates the coefficients for a filter that has a frequency response given by a
Butterworth filter with order `N` and cutoff frequency `freq_c`, squared. 

Uses `N` exponentials and `N/2` sets of coefficients (a,b,c,d). N should therefore be even,
since exponentials come in pairs to ensure a real-valued filter.

However, the special case N=1 is allowed, which gives a single (real) exponential filter.

Frequency response: `Ghat(omega) = 1 / (1 + (omega / freq_c)^(2*N))`
Real filter shape: `G(t) = sum_{i=1}^{N/2} exp(-c_i*abs(t))*(a_i*cos(d_i * abs(t)) + b_i*sin(d_i * abs(t)))`

This function supports two types of filters:
* A **single exponential filter** when `N=1`. This is a special case that
  generates two coefficients instead of 4. The unidirectional filter is a single exponential,
  and `N_coeffs = 0.5`. Only `a1` and `c1` are returned.
* A **Butterworth squared filter** for `N>1`. This generates `N/2` sets of
  coefficients (a,b,c,d), representing a filter of order `N`. The
  coefficients are computed based on the filter's order and cutoff frequency.

Arguments
=========
- `N`: The order parameter for the filter. `N=1` for a single exponential.
  For `N>1`, the filter's order is `N`. Must be a non-negative even integer.
- `freq_c`: The cutoff frequency of the filter. Must be a real number.

Returns
=======
- A `NamedTuple` containing the filter coefficients and `N_coeffs`, the number
  of coefficient pairs.
"""
function set_offline_BW2_filter_params(;N::Int=1,freq_c::Real=1) 
    if N == 1
        # Special case single exponential
        N_coeffs = 0.5
    elseif N/2 != floor(N/2) || N < 0 
        error("N must be a non-negative even integer, or 1 for a single exponential.")
    else
        N_coeffs = Int(N/2)
    end

    freq_c = abs(freq_c) # Ensure freq_c is positive
    if N_coeffs == 0.5 # special case N=0, single exponential only has a cosine component
        a1 = freq_c/2 
        c1 = freq_c
        filter_params = (; a1 = a1, c1 = c1, N_coeffs = N_coeffs)
        return filter_params

    else
        filter_params = NamedTuple()
        for i in 1:N_coeffs
            
            a = (freq_c/N)*sin(pi/2/N*(2*i-1))
            b = (freq_c/N)*cos(pi/2/N*(2*i-1))
            c = freq_c*sin(pi/2/N*(2*i-1))
            d = freq_c*cos(pi/2/N*(2*i-1))

            temp_params = NamedTuple{(Symbol("a$i"), Symbol("b$i"),Symbol("c$i"),Symbol("d$i"))}([a,b,c,d])
            filter_params = merge(filter_params,temp_params)
        end

        return merge(filter_params, (; N_coeffs = N_coeffs))
    end
end

"""
    set_online_BW_filter_params(; N::Int=1, freq_c::Real=1)

Calculates the coefficients for a filter that has a frequency response given by a
Butterworth filter with order `N` and cutoff frequency `freq_c`. Note that the frequency
response is not squared, like in the offline forward-backward filter, and the frequency response
is not real-valued, implying a nonlinear phase shift. 

Uses `N` exponentials and `N/2` sets of coefficients (a,b,c,d). N should therefore be even,
since exponentials come in pairs to ensure a real-valued filter.

However, the special case N=1 is allowed, which gives a single (real) exponential filter.

Frequency response: `abs(Ghat(omega)) = 1 / sqrt(1 + (omega / freq_c)^(2*N))`
Real filter shape: `G(t) = sum_{i=1}^{N/2} exp(-c_i*t)* (a_i*cos(d_i * t) + b_i*sin(d_i * t))` 
for t>=0, and 0 for t<0.

This function supports two types of filters:
* A **single exponential filter** when `N=1`. This is a special case that
  generates two coefficients instead of 4. The unidirectional filter is a single exponential,
  and `N_coeffs = 0.5`. Only `a1` and `c1` are returned.
* A **Butterworth filter** for `N>1`. This generates `N/2` sets of
  coefficients (a,b,c,d), representing a filter of order `N`. The
  coefficients are computed based on the filter's order and cutoff frequency.

Arguments
=========
- `N`: The order parameter for the filter. `N=1` for a single exponential.
  For `N>1`, the filter's order is `N`. Must be a non-negative even integer.
- `freq_c`: The cutoff frequency of the filter. Must be a real number.

Returns
=======
- A `NamedTuple` containing the filter coefficients and `N_coeffs`, the number
  of coefficient pairs.
"""
function set_online_BW_filter_params(;N::Int=1,freq_c::Real=1) 
    if N == 1
        # Special case single exponential
        N_coeffs = 0.5
    elseif N/2 != floor(N/2) || N < 0 
        error("N must be a non-negative even integer, or 1 for a single exponential.")
    else
        N_coeffs = Int(N/2)
    end

    freq_c = abs(freq_c) # Ensure freq_c is positive
    if N_coeffs == 0.5 # special case N=0, single exponential only has a cosine component
        a1 = freq_c
        c1 = freq_c
        filter_params = (; a1 = a1, c1 = c1, N_coeffs = N_coeffs)
        return filter_params

    else
        filter_params = NamedTuple()
        for i in 1:N_coeffs
            
            c = freq_c*sin(pi/2/N*(2*i-1))
            d = -freq_c*cos(pi/2/N*(2*i-1))
            ri = exp(pi*im*(2*i-1)/ (2*N)) # ith of the 2N roots of -1
            Ai = 1
            for k in 1:N
                if k != i
                    rk = exp(pi*im*(2*k-1)/ (2*N))
                    Ai *= 1/(ri - rk)
                end
            end
            b = 2*freq_c*real(Ai*exp(-im*pi*N/2))
            a = -2*freq_c*imag(Ai*exp(-im*pi*N/2))

            temp_params = NamedTuple{(Symbol("a$i"), Symbol("b$i"),Symbol("c$i"),Symbol("d$i"))}([a,b,c,d])
            filter_params = merge(filter_params,temp_params)
        end

        return merge(filter_params, (; N_coeffs = N_coeffs))
    end
end

"""
    create_original_vars(config::AbstractConfig)

Creates a `NamedTuple` to serve as auxiliary fields for the original variables 
in a simulation. The fields are initialised from the input data to ensure that
they are in the correct location. - the value of the data is unimportant.

Arguments
=========
- `config`: An instance of `AbstractConfig` containing the names of the
  variables and the simulation grid.

Returns
=======
A `NamedTuple` where each key is a `Symbol` of a variable name to be filtered,
and each value is an empty `CenterField` for that variable.
"""
function create_original_vars(config::AbstractConfig)

    var_names_to_filter = config.var_names_to_filter
    grid = config.grid
    architecture = config.architecture
    backend = config.backend
    vars = Dict()
    original_data_filename = config.original_data_filename

    for var_name in var_names_to_filter
        fts_data = FieldTimeSeries(original_data_filename, var_name, architecture=architecture,backend=backend)[1]
        vars[Symbol(var_name)] = fts_data
    end
    return NamedTuple(vars)
end

"""
    create_filtered_vars(config::AbstractConfig)

Creates a `Tuple` of `Symbol`s representing the names of the filtered tracer
variables.

* For a single-exponential filter (`N_coeffs = 0.5`), the function generates
  names with a `_C1` suffix.
* For a multi-coefficient filter (`N_coeffs > 0.5`), it generates pairs of
  names for each coefficient, suffixed with `_C#` and `_S#`, where `#` is the
  coefficient index.

If `map_to_mean` or `compute_mean_velocities` is enabled in the configuration,
additional symbols are created for the spatial mapping variables corresponding
to each velocity component, prefixed with `xi_` and suffixed with the corresponding
coefficient names.

Arguments
=========
- `config`: An instance of `AbstractConfig` containing the names of the variables
  to filter, the filter parameters, and the `map_to_mean` and `compute_mean_velocities` booleans.

Returns
=======
A `Tuple` of `Symbol`s representing the names of the filtered variables to be
used as tracers in the simulation.
"""
function create_filtered_vars(config::AbstractConfig)
    
    var_names_to_filter = config.var_names_to_filter
    velocity_names = config.velocity_names
    filter_params = config.filter_params
    map_to_mean = config.map_to_mean
    compute_mean_velocities = config.compute_mean_velocities
    N_coeffs = filter_params.N_coeffs
    label = config.label

    if N_coeffs == 0.5 # special case, single exponential only has a cosine component
        gC_symbols = Symbol[]
        for var_name in var_names_to_filter
            push!(gC_symbols, Symbol(var_name, label, "_C1"))
        end
        # May also need xi maps. We need one for every velocity dimension, so lets use the velocity names to name them
        if map_to_mean || compute_mean_velocities
            for vel_name in velocity_names
                push!(gC_symbols, Symbol("xi_", vel_name, label,"_C1"))
            end
        end

        return Tuple(gC_symbols)

    else
        gC_symbols = Symbol[]
        gS_symbols = Symbol[]

        for var_name in var_names_to_filter
            for i in 1:N_coeffs
                push!(gC_symbols, Symbol(var_name, label, "_C", i))
                push!(gS_symbols, Symbol(var_name, label, "_S", i))
            end
        end

        # May also need xi maps. We need one for every velocity dimension, so lets use the velocity names to name them
        if map_to_mean || compute_mean_velocities
            for vel_name in velocity_names
                for i in 1:N_coeffs
                    push!(gC_symbols, Symbol("xi_", vel_name, label, "_C", i))
                    push!(gS_symbols, Symbol("xi_", vel_name, label, "_S", i))
                end
            end
        end

        return (Tuple(gC_symbols)..., Tuple(gS_symbols)...)
    end
end

"""
    _make_gC_forcing(i::Int, var_name::String, filter_params::NamedTuple)

Create a forcing term for the cosine component (gC) of a filtered variable.
Includes the special case of a single exponential.

# Arguments
- `i::Int`: The index of the coefficient pair (cᵢ, dᵢ) to use from `filter_params`.
- `labelled_var_name::String`: The name of the variable being filtered (e.g., "T")
    including label if used.
- `filter_params::NamedTuple`: A `NamedTuple` containing all filter coefficients.

# Returns
- A `Forcing` object configured to compute the forcing term for the gC field.
"""
function _make_gC_forcing(i::Int, labelled_var_name::String, filter_params::NamedTuple)
    if filter_params.N_coeffs == 0.5 # Single exponential special case has a simpler forcing
        c = getproperty(filter_params, Symbol("c",i))
        gCkey = Symbol(labelled_var_name,"_C",i)
        forcing_func = (args...) -> -args[end][1]*args[end-1] 
        return Forcing(forcing_func, parameters = (c,), field_dependencies = (gCkey,))
    else
        c = getproperty(filter_params, Symbol("c",i))
        d = getproperty(filter_params, Symbol("d",i))
        gCkey = Symbol(labelled_var_name, "_C",i)
        gSkey = Symbol(labelled_var_name, "_S",i)
        forcing_func = (args...) -> -args[end][1]*args[end-2] - args[end][2]*args[end-1]
        return Forcing(forcing_func, parameters = (c,d), field_dependencies = (gCkey,gSkey))
    end
end

"""
    _make_gS_forcing(i::Int, var_name::String, filter_params::NamedTuple)

Create a forcing term for the sine component (gS) of a filtered variable.

# Arguments
- `i::Int`: The index of the coefficient pair (cᵢ, dᵢ) to use from `filter_params`.
- `labelled_var_name::String`: The name of the variable being filtered (e.g., "T")
    including label if used.
- `filter_params::NamedTuple`: A `NamedTuple` containing all filter coefficients.

# Returns
- A `Forcing` object configured to compute the forcing term for the gS field.
"""
function _make_gS_forcing(i::Int, labelled_var_name::String, filter_params::NamedTuple)
    c = getproperty(filter_params, Symbol("c",i))
    d = getproperty(filter_params, Symbol("d",i))
    gCkey = Symbol(labelled_var_name,"_C",i)
    gSkey = Symbol(labelled_var_name,"_S",i)  
    forcing_func = (args...) -> -args[end][1]*args[end-1] + args[end][2]*args[end-2]
    return Forcing(forcing_func, parameters = (c,d), field_dependencies = (gCkey,gSkey))
end

"""
    _make_xiC_forcing(i::Int, vel_name::String, filter_params::NamedTuple)

Create a forcing term for the cosine component (xiC) of a map variable. The function handles 
the special case of a single exponential filter where `d` is zero.

# Arguments
- `i::Int`: The index of the coefficient pair (cᵢ, dᵢ) to use from `filter_params`.
- `vel_name::String`: The name of the velocity variable (e.g., "u").
- `filter_params::NamedTuple`: A `NamedTuple` containing all filter coefficients.

# Returns
- A `Forcing` object configured to compute the forcing term for the xiC field.
"""
function _make_xiC_forcing(i::Int, vel_name::String, filter_params::NamedTuple)
    c = getproperty(filter_params, Symbol("c",i))
    d = filter_params.N_coeffs == 0.5 ? 0 : getproperty(filter_params, Symbol("d",i)) # Single exponential special case gets d = 0
    forcing_func = (args...) -> -args[end][1]/(args[end][1]^2 + args[end][2]^2)*args[end-1]
    return Forcing(forcing_func, parameters = (c,d), field_dependencies = (Symbol(vel_name),))
end

"""
    _make_xiS_forcing(i::Int, vel_name::String, filter_params::NamedTuple)

Create a forcing term for the sine component (xiS) of a map variable. 

# Arguments
- `i::Int`: The index of the coefficient pair (cᵢ, dᵢ) to use from `filter_params`.
- `vel_name::String`: The name of the velocity variable (e.g., "u").
- `filter_params::NamedTuple`: A `NamedTuple` containing all filter coefficients.

# Returns
- A `Forcing` object configured to compute the forcing term for the xiS field.
"""
function _make_xiS_forcing(i::Int, vel_name::String, filter_params::NamedTuple)
    c = getproperty(filter_params, Symbol("c",i))
    d = getproperty(filter_params, Symbol("d",i))
    forcing_func = (args...) -> -args[end][2]/(args[end][1]^2 + args[end][2]^2)*args[end-1]
    return Forcing(forcing_func, parameters = (c,d), field_dependencies = (Symbol(vel_name),))
end


"""
    create_forcing(filtered_vars::Tuple{Vararg{Symbol}}, config::AbstractConfig)

Creates a `NamedTuple` of forcing functions for each filtered variable and,
if enabled, for the spatial mapping variables. These forcing terms are used
to numerically integrate the filter equations.

The function handles two cases: a single-exponential filter
(`N_coeffs = 0.5`) and a multi-coefficient Butterworth squared filter
(`N_coeffs > 0.5`).

* For standard filtered variables, the forcing is a combination of terms
  derived from the filter's coefficients and a term from the original data.
* For spatial mapping variables (if `map_to_mean` or `compute_mean_velocities` is true), the forcing
  includes terms derived from the filter's coefficients and a term from the
  original velocity data.

Arguments
=========
- `filtered_vars`: A `Tuple` of `Symbol`s representing the names of the
  filtered variables.
- `config`: An instance of `AbstractConfig` containing the names of the
  variables to be filtered, velocity names, and the filter parameters.

Returns
=======
A `NamedTuple` where each key is a variable name from `filtered_vars` and
each value is a `Tuple` of the corresponding forcing functions.
"""
function create_forcing(filtered_vars::Tuple{Vararg{Symbol}}, config::AbstractConfig)

    var_names_to_filter = config.var_names_to_filter
    velocity_names = config.velocity_names
    filter_params = config.filter_params
    N_coeffs = filter_params.N_coeffs
    label = config.label
    # Initialize dictionary
    gC_forcings_dict = Dict()
    

    # Make a simple forcing function for original data forcing - the final argument is the field dependence.
    original_var_forcing_func(args...) = args[end]

    # Special case `N_coeffs = 0.5`, single exponential
    if N_coeffs == 0.5
        # Build by hand as only one coefficient
        for var_name in var_names_to_filter
            labelled_var_name = var_name * label
            var_key = Symbol(var_name)
            gCkey = Symbol(labelled_var_name,"_C1")

            # The forcing for gC is the sum of a filter forcing term and the original data forcing
            gC_forcing = _make_gC_forcing(1, labelled_var_name, filter_params)
            gC_original_var_forcing = Forcing(original_var_forcing_func,field_dependencies = (;var_key))
            gC_forcings_dict[gCkey] = (gC_forcing, gC_original_var_forcing)
        end

        # Check if we need xi forcing (implied by number of filtered_vars)
        if length(filtered_vars) > N_coeffs*length(var_names_to_filter)*2

            # Create forcing for xi maps (one for each velocity component)
            for vel_name in velocity_names
                labelled_var_name = "xi_" * vel_name * label
                gCkey = Symbol(labelled_var_name, "_C1") 
                 
                # The forcing for xiC includes a term involving xiC (as for tracers) and a 
                # term involving the corresponding velocity
                gC_forcing = _make_gC_forcing(1, labelled_var_name, filter_params)
                xiC_forcing = _make_xiC_forcing(1, vel_name, filter_params)
                gC_forcings_dict[gCkey] = (xiC_forcing, gC_forcing)

            end
        end

        return NamedTuple(gC_forcings_dict)

    else
        gS_forcings_dict = Dict()

        for var_name in var_names_to_filter
            var_key = Symbol(var_name)
            for i in 1:N_coeffs
                labelled_var_name = var_name * label
                gCkey = Symbol(labelled_var_name,"_C",i)
                gSkey = Symbol(labelled_var_name,"_S",i)

                # The forcing for gC is the sum of a filter forcing term and the original data forcing
                gC_forcing_i = _make_gC_forcing(i, labelled_var_name, filter_params)
                gC_original_var_forcing = Forcing(original_var_forcing_func, field_dependencies= (;var_key))
                gC_forcings_dict[gCkey] = (gC_forcing_i, gC_original_var_forcing)

                # The forcing for gS is just a filter forcing term
                gS_forcings_dict[gSkey] = _make_gS_forcing(i, labelled_var_name, filter_params)

            end
        end

        # Check if we need xi forcing (implied by number of filtered_vars)
        if length(filtered_vars) > N_coeffs*length(var_names_to_filter)*2
            # Create forcing for xi maps (one for each velocity component)
            for vel_name in velocity_names
                
                for i in 1:N_coeffs
                    labelled_var_name = "xi_" * vel_name * label
                    gCkey = Symbol(labelled_var_name, "_C", i) 
                    gSkey = Symbol(labelled_var_name, "_S", i) 

                    # The forcing for xiC includes a term involving xiC and xiS (as for tracers, so reuse forcing constructor) and a 
                    # term involving the corresponding velocity
                    gC_forcing_i = _make_gC_forcing(i, labelled_var_name, filter_params)
                    xiC_forcing_i = _make_xiC_forcing(i, vel_name, filter_params)
                    gC_forcings_dict[gCkey] = (xiC_forcing_i, gC_forcing_i)

                    # The forcing for xiS also includes a term involving xiC and xiS and a term involving the corresponding velocity
                    gS_forcing_i = _make_gS_forcing(i, labelled_var_name, filter_params)
                    xiS_forcing_i = _make_xiS_forcing(i, vel_name, filter_params)
                    gS_forcings_dict[gSkey] = (xiS_forcing_i, gS_forcing_i)

                end
            end
        end

        return (; NamedTuple(gC_forcings_dict)..., NamedTuple(gS_forcings_dict)...)
    end
end

"""
    create_output_fields(model::AbstractModel, config::AbstractConfig)

Reconstructs the final output fields from the model's tracers and auxiliary
fields. This function performs the following steps:

1.  **Reconstructs filtered variables**: For each variable to be filtered, it sums
    the contributions from the individual filter coefficients (`_C` and `_S`
    tracers) using the coefficients from `filter_params`.
2.  **Reconstructs spatial mapping fields**: If spatial mapping is enabled,
    the function also reconstructs the `xi_` fields that represent the filtered
    position.
3.  **Builds mean velocities**: If `compute_mean_velocities` is true, the function
    reconstructs the mean velocity fields using the `xi_` fields.
4.  **Includes original data**: The original data is added to the output
    dictionary for comparison and analysis if `config.output_original_data`
    is true.

Arguments
=========
- `model`: An instance of an `AbstractModel` containing the tracer and
  auxiliary fields.
- `config`: An instance of `AbstractConfig` with the names of the variables,
  velocity components, and filter parameters.

Returns
=======
A `Dict` where keys are the names of the output fields (e.g.,
`var_name_Lagrangian_filtered`, `xi_vel_name`, `var_name`) and values are the
corresponding reconstructed `Field`s.
"""
function create_output_fields(model::AbstractModel, config::AbstractConfig)

    var_names_to_filter = config.var_names_to_filter
    velocity_names = config.velocity_names
    filter_params = config.filter_params
    N_coeffs = filter_params.N_coeffs
    map_to_mean = config.map_to_mean
    compute_mean_velocities = config.compute_mean_velocities
    label = config.label
    outputs_dict = Dict()

    # When offline filtering, we can turn off advection to get Eulerian filtered fields
    if config.filter_mode == "offline" && config.advection === nothing
        filter_identifier = "_Eulerian_filtered"
    else
        filter_identifier = "_Lagrangian_filtered"
    end

    for var_name in var_names_to_filter
        labelled_var_name = var_name * label
        if N_coeffs == 0.5
            # Special case, single exponential only has a cosine component
            gC1 = getproperty(model.tracers,Symbol(labelled_var_name * "_C1"))
            g_total = filter_params.a1 * gC1
            outputs_dict[labelled_var_name * filter_identifier] = g_total
        else
            # Reconstruct the filtered tracer fields, starting with the first coefficient
            gC1 = getproperty(model.tracers, Symbol(labelled_var_name * "_C1"))
            gS1 = getproperty(model.tracers, Symbol(labelled_var_name * "_S1"))
            g_total = filter_params.a1 * gC1 + filter_params.b1 * gS1

            # Then add the other coefficients
            for i in 2:N_coeffs
                a = getproperty(filter_params,Symbol("a$i"))
                b = getproperty(filter_params,Symbol("b$i"))
                gCi = getproperty(model.tracers,Symbol(labelled_var_name * "_C$i" ))
                gSi = getproperty(model.tracers,Symbol(labelled_var_name * "_S$i" ))
                g_total += a * gCi + b * gSi
            end
            outputs_dict[labelled_var_name * filter_identifier] = g_total
        end
    end

    # Reconstruct the maps, if we map to mean
    if map_to_mean
        for vel_name in velocity_names
            labelled_var_name = "xi_" * vel_name * label
            if N_coeffs == 0.5
                # Special case, single exponential only has a cosine component
                xiC1 = getproperty(model.tracers,Symbol(labelled_var_name * "_C1"))
                g_total = filter_params.a1 * xiC1
                outputs_dict[labelled_var_name] = g_total
            else
                # Start with the first coefficient
                xiC1 = getproperty(model.tracers,Symbol(labelled_var_name * "_C1"))
                xiS1 = getproperty(model.tracers,Symbol(labelled_var_name * "_S1"))
                g_total = filter_params.a1 * xiC1 + filter_params.b1 * xiS1

                # Then add the other coefficients
                for i in 2:N_coeffs
                    a = getproperty(filter_params,Symbol("a$i"))
                    b = getproperty(filter_params,Symbol("b$i"))
                    xiCi = getproperty(model.tracers,Symbol(labelled_var_name * "_C$i"))
                    xiSi = getproperty(model.tracers,Symbol(labelled_var_name * "_S$i"))
                    g_total += a * xiCi + b * xiSi
                end
                outputs_dict[labelled_var_name] = g_total
            end
        end
    end

    # Reconstruct the mean velocities
    if compute_mean_velocities
        for vel_name in velocity_names
            labelled_var_name = "xi_" * vel_name * label
            if N_coeffs == 0.5
                # Special case, single exponential only has a cosine component
                xiC1 = getproperty(model.tracers,Symbol(labelled_var_name * "_C1"))
                g_total = - filter_params.a1 * filter_params.c1 * xiC1
                outputs_dict[vel_name * label * filter_identifier] = g_total
            else
                # Start with the first coefficient
                xiC1 = getproperty(model.tracers,Symbol(labelled_var_name * "_C1"))
                xiS1 = getproperty(model.tracers,Symbol(labelled_var_name * "_S1"))
                g_total = ((-filter_params.a1 *filter_params.c1 + filter_params.b1 * filter_params.d1) * xiC1 
                + (-filter_params.a1 * filter_params.d1 - filter_params.b1 * filter_params.c1) * xiS1)

                # Then add the other coefficients
                for i in 2:N_coeffs
                    a = getproperty(filter_params,Symbol("a$i"))
                    b = getproperty(filter_params,Symbol("b$i"))
                    c = getproperty(filter_params,Symbol("c$i"))
                    d = getproperty(filter_params,Symbol("d$i"))
                    xiCi = getproperty(model.tracers,Symbol(labelled_var_name * "_C$i"))
                    xiSi = getproperty(model.tracers,Symbol(labelled_var_name * "_S$i"))
                    g_total += (-a * c + b * d) * xiCi + (-a * d - b * c) * xiSi
                end
                outputs_dict[vel_name * label * filter_identifier] = g_total
            end
        end
    end

    # We can also add the saved vars for comparison if this is an offline filter, otherwise do this manually 
    if config.filter_mode == "offline" && config.output_original_data
        for var_name in var_names_to_filter
            outputs_dict[var_name] = getproperty(model.auxiliary_fields, Symbol(var_name))
        end
        for vel_name in velocity_names
            outputs_dict[vel_name] = getproperty(model.velocities, Symbol(vel_name))
        end
    end
 

    return outputs_dict
end

"""
    update_input_data!(model::AbstractModel, input_data::NamedTuple)

Updates the velocity and auxiliary fields of a simulation at the current
simulation time `t`. This function is designed to be used as a callback in an
Oceananigans `Simulation` at callsite UpdateStateCallsite().

The function performs two main tasks:
1.  **Updates velocities**: It sets the `u`, `v`, and `w` velocity fields of
    the `model` to the corresponding data from the `velocity_data`
    `FieldTimeSeries` at the current simulation time.
2.  **Updates auxiliary fields**: It updates the auxiliary fields of the
    `model` with the original data from the `var_data` `FieldTimeSeries`, which
    are used for forcing terms.

Arguments
=========

- `model`: The model.
- `input_data`: A `NamedTuple` containing `velocity_data` and `var_data`,
  where each field is a `Tuple` of `FieldTimeSeries` objects.
"""
function update_input_data!(model::AbstractModel, input_data::NamedTuple)
    velocity_timeseries = input_data.velocity_data
    original_var_timeseries = input_data.var_data
    t = model.clock.time
    
    # Update the velocities
    kwargs = (; (Symbol(vel_fts.name) => vel_fts[Time(t)] for vel_fts in velocity_timeseries)...)
    set!(model; kwargs...)          

    # We also update the saved original variables to be used for forcing - these are auxiliary fields so need to be set separately
    for original_var_fts in original_var_timeseries
        set!(getproperty(model.auxiliary_fields, Symbol(original_var_fts.name)), original_var_fts[Time(t)])
        # halo regions get filled automatically
    end
end

"""
    initialise_filtered_vars_from_data(model::AbstractModel, saved_original_vars::Tuple,
                             config::AbstractConfig)

Initializes the model's tracer fields, which represent the components of the
filtered variables. This function sets the initial values of the filtered
variables to the (scaled) first timestep of the original data. This improves
the "spin-up" of the filter simulation by providing a good starting point.

The initialization formula depends on the number of filter coefficients
(`N_coeffs`):

- For a **single-exponential filter** (`N_coeffs = 0.5`), only the `_C1`
  tracer exists and is initialized.
- For a **multi-coefficient filter** (`N_coeffs > 0.5`), both the `_C` and `_S`
  tracers for each coefficient are initialized.

Arguments
=========
- `model`: The `AbstractModel` whose tracers are to be initialized.
- `saved_original_vars`: A `Tuple` of `FieldTimeSeries` objects containing
  the original data for each variable.
- `config`: An instance of `AbstractConfig` with the filter parameters.
"""
function initialise_filtered_vars_from_data(model::AbstractModel, input_data::NamedTuple, config::AbstractConfig)
    filter_params = config.filter_params
    label = config.label
    for original_var_fts in input_data.var_data
        var_name = original_var_fts.name
        labelled_var_name = var_name * label
        if filter_params.N_coeffs == 0.5 # Special case of single exponential
            filtered_var_C = Symbol(labelled_var_name,"_C1",)
            c1 = filter_params.c1
            set!(getproperty(model.tracers, filtered_var_C), 1/c1*original_var_fts[Time(0)])
        else
            for i in 1:filter_params.N_coeffs
                filtered_var_C = Symbol(labelled_var_name,"_C",i)
                filtered_var_S = Symbol(labelled_var_name,"_S",i)
                ci = getproperty(filter_params,Symbol("c$i"))
                di = getproperty(filter_params,Symbol("d$i"))
                set!(getproperty(model.tracers, filtered_var_C), ci/(ci^2 + di^2)*original_var_fts[Time(0)])
                set!(getproperty(model.tracers, filtered_var_S), di/(ci^2 + di^2)*original_var_fts[Time(0)])
            end
        end
    end
end

"""
    initialise_filtered_vars_from_model(model::AbstractModel,config::AbstractConfig)


Initializes the model's filtered tracer fields using the actual tracer fields that are
assumed to have been already set. This improves the "spin-up" of the filter simulation 
by providing a good starting point.

The initialization formula depends on the number of filter coefficients
(`N_coeffs`):

- For a **single-exponential filter** (`N_coeffs = 0.5`), only the `_C1`
  tracer exists and is initialized.
- For a **multi-coefficient filter** (`N_coeffs > 0.5`), both the `_C` and `_S`
  tracers for each coefficient are initialized.

Arguments
=========
- `model`: The `AbstractModel` whose tracers are to be initialized.
- `config`: An instance of `AbstractConfig` with the filter parameters.
"""
function initialise_filtered_vars_from_model(model::AbstractModel, config::AbstractConfig)
    filter_params = config.filter_params
    var_names_to_filter = config.var_names_to_filter
    label = config.label
    for var_name in var_names_to_filter
        labelled_var_name = var_name * label
        if filter_params.N_coeffs == 0.5 # Special case of single exponential
            filtered_var_C = Symbol(labelled_var_name,"_C1",)
            c1 = filter_params.c1
            # original data can be tracer or auxiliary field
            if Symbol(var_name) in propertynames(model.tracers)
                field = getproperty(model.tracers,Symbol(var_name))
            elseif Symbol(var_name) in propertynames(model.auxiliary_fields)
                field = getproperty(model.auxiliary_fields,Symbol(var_name))
            else
                error("Variable $var_name not found in model tracers or auxiliary fields.")
            end
            set!(getproperty(model.tracers, filtered_var_C), 1/c1*field)
        else
            for i in 1:filter_params.N_coeffs
                filtered_var_C = Symbol(labelled_var_name,"_C",i)
                filtered_var_S = Symbol(labelled_var_name,"_S",i)
                ci = getproperty(filter_params,Symbol("c$i"))
                di = getproperty(filter_params,Symbol("d$i"))
                if Symbol(var_name) in propertynames(model.tracers)
                    field = getproperty(model.tracers,Symbol(var_name))
                elseif Symbol(var_name) in propertynames(model.auxiliary_fields)
                    field = getproperty(model.auxiliary_fields,Symbol(var_name))
                else
                    error("Variable $var_name not found in model tracers or auxiliary fields.")
                end
                set!(getproperty(model.tracers, filtered_var_C), ci/(ci^2 + di^2)*field)
                set!(getproperty(model.tracers, filtered_var_S), di/(ci^2 + di^2)*field)
            end
        end
    end
end

"""
    zero_closure_for_filtered_vars(config::AbstractConfig)


Initializes the model's filtered tracer fields using the actual tracer fields that are
assumed to have been already set. This improves the "spin-up" of the filter simulation 
by providing a good starting point.

The initialization formula depends on the number of filter coefficients
(`N_coeffs`):

- For a **single-exponential filter** (`N_coeffs = 0.5`), only the `_C1`
  tracer exists and is initialized.
- For a **multi-coefficient filter** (`N_coeffs > 0.5`), both the `_C` and `_S`
  tracers for each coefficient are initialized.

Arguments
=========
- `model`: The `AbstractModel` whose tracers are to be initialized.
- `config`: An instance of `AbstractConfig` with the filter parameters.
"""
function zero_closure_for_filtered_vars(config::AbstractConfig)
    var_names_to_filter = config.var_names_to_filter
    N_coeffs = config.filter_params.N_coeffs
    map_to_mean = config.map_to_mean
    compute_mean_velocities = config.compute_mean_velocities
    label = config.label
    dict = Dict()
    for var_name in var_names_to_filter
        labelled_var_name = var_name * label
        if N_coeffs == 0.5 # Special case of single exponential
            filtered_var_C = Symbol(labelled_var_name,"_C1",)
            dict[filtered_var_C] = 0.0
        else
            for i in 1:N_coeffs
                filtered_var_C = Symbol(labelled_var_name,"_C",i)
                filtered_var_S = Symbol(labelled_var_name,"_S",i)
                dict[filtered_var_C] = 0.0
                dict[filtered_var_S] = 0.0
            end
        end
    end
    if map_to_mean || compute_mean_velocities
        velocity_names = config.velocity_names
        for vel_name in velocity_names
            if N_coeffs == 0.5 # Special case of single exponential
                filtered_var_C = Symbol("xi_", vel_name, label, "_C1",)
                dict[filtered_var_C] = 0.0
            else
                for i in 1:N_coeffs
                    filtered_var_C = Symbol("xi_", vel_name, label, "_C",i)
                    filtered_var_S = Symbol("xi_", vel_name, label, "_S",i)
                    dict[filtered_var_C] = 0.0
                    dict[filtered_var_S] = 0.0
                end
            end
        end
    end
    filtered_closure = NamedTuple(dict)
    return filtered_closure
end


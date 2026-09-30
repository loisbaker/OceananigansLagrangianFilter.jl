using Oceananigans: AbstractModel

"""
    copy_file_metadata!(original_file::JLD2.JLDFile, new_file::JLD2.JLDFile,
                        timeseries_vars_to_copy::Tuple{Vararg{String}})

Copies essential metadata from an existing Oceananigans JLD2 output file to a
new file.

This function is a utility for creating a new file structure with all
necessary simulation metadata (such as grid information and serialized
objects) without copying large timeseries data. It copies core simulation
metadata (`serialized`) and the `serialized` entries for
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
    _copy_jld2_recursive!(original_file, new_file, "serialized")
    

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
- A `NamedTuple` containing the filter coefficients, `N_coeffs` (the number
  of coefficient pairs), and the output options `outputs = :combined` and
  `sine_parity = :even` (the defaults, see `filter_output_options`).
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
    if N_coeffs == 0.5 # special case N=1, single exponential only has a cosine component
        a1 = freq_c/2 
        c1 = freq_c
        filter_params = (; a1 = a1, c1 = c1, N_coeffs = N_coeffs, outputs = :combined, sine_parity = :even)
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

        return merge(filter_params, (; N_coeffs = N_coeffs, outputs = :combined, sine_parity = :even))
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
- A `NamedTuple` containing the filter coefficients, `N_coeffs` (the number
  of coefficient pairs), and the output option `outputs = :combined` (the default,
  see `filter_output_options`).
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
    if N_coeffs == 0.5 # special case N=1, single exponential only has a cosine component
        a1 = freq_c
        c1 = freq_c
        filter_params = (; a1 = a1, c1 = c1, N_coeffs = N_coeffs, outputs = :combined)
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

        return merge(filter_params, (; N_coeffs = N_coeffs, outputs = :combined))
    end
end

"""
    set_offline_spectrum_filter_params(; alpha::Real, freqs::AbstractVector, normalisation::Symbol = :spectral)

Coefficients for an exponentially-windowed spectral filter, which extracts the signal in a narrow band
around each frequency in `freqs`. The window is `w(t) = exp(-alpha*|t|)`, giving the in-phase (`C`) and
quadrature (`S`) weight functions at each frequency `omega_n`:

    C(t) = A_n exp(-alpha*|t|) cos(omega_n*t)
    S(t) = A_n exp(-alpha*|t|) sin(omega_n*t)

so that `a_n = b_n = A_n`, `c_n = alpha` and `d_n = omega_n`. The `S` weight functions are odd in `t`
(`sine_parity = :odd`), and each term is output separately (`outputs = :separate`).

`alpha` sets the frequency resolution: the window lasts about `1/alpha`, so frequencies closer than
about `alpha` are not resolved.

Two normalisations are available:
- `:spectral` (default): `A_n = sqrt(alpha)`, so that the window has unit energy, `∫ (A_n w)^2 dt = 1`.
  Then `C^2 + S^2` estimates the (two-sided) power spectral density of the signal at `omega_n`,
  normalised so that the variance is `∫ S(omega) domega / 2π`.
- `:unit_gain`: `A_n = alpha (alpha^2 + 4omega_n^2) / (2alpha^2 + 4omega_n^2)`, so that `C` passes a signal at
  exactly `omega_n` through unchanged: `C` is the band-passed signal, and `S` is its quadrature (a quarter
  period later). Then `sqrt(C^2 + S^2)` is its envelope and `atan(S, C)` its instantaneous phase (`omega_n*t`
  plus the phase of the signal). `S` has gain `4omega_n^2 / (2alpha^2 + 4omega_n^2)`, so the envelope and phase
  are accurate to about `alpha^2 / (2omega_n^2)`. At `omega_n = 0` this reduces to `alpha/2`, the single
  exponential of `set_offline_BW2_filter_params(N = 1)`.

The `S` terms can be left out of the outputs by setting their coefficients `b_n` to zero. See
[`set_online_spectrum_filter_params`](@ref) for the online filter.

Arguments
=========
- `alpha`: Decay rate of the exponential window. Must be positive.
- `freqs`: Frequencies (radians per unit time) at which to extract the signal.
- `normalisation`: `:spectral` (default) or `:unit_gain`.

Returns
=======
- A `NamedTuple` of coefficients, `N_coeffs`, `outputs = :separate` and `sine_parity = :odd`.
"""
function set_offline_spectrum_filter_params(; alpha::Real, freqs::AbstractVector, normalisation::Symbol = :spectral)
    amplitude(omega) = normalisation === :spectral ? sqrt(alpha) :
                       alpha * (alpha^2 + 4omega^2) / (2alpha^2 + 4omega^2) # Unit gain at omega
    filter_params = spectrum_filter_coefficients(amplitude; alpha, freqs, normalisation)
    return merge(filter_params, (; sine_parity = :odd))
end

"""
    set_online_spectrum_filter_params(; alpha::Real, freqs::AbstractVector, normalisation::Symbol = :spectral)

Coefficients for an exponentially-windowed spectral filter for the online filter: the causal version of
[`set_offline_spectrum_filter_params`](@ref), which only uses the past. The window is `w(t) = exp(-alpha*t)`
for `t > 0` (the time before the present), giving the in-phase (`C`) and quadrature (`S`) weight functions
at each frequency `omega_n`:

    C(t) = A_n exp(-alpha*t) cos(omega_n*t)
    S(t) = A_n exp(-alpha*t) sin(omega_n*t)

so that `a_n = b_n = A_n`, `c_n = alpha` and `d_n = omega_n`, and each term is output separately
(`outputs = :separate`).

`alpha` sets the frequency resolution: the window lasts about `1/alpha`, so frequencies closer than
about `alpha` are not resolved. Because the window only uses the past, the envelope `sqrt(C^2 + S^2)`
responds to changes in the signal with a delay of about `1/alpha`.

Two normalisations are available:
- `:spectral` (default): `A_n = sqrt(2alpha)`, so that the window has unit energy, `∫ (A_n w)^2 dt = 1`.
  Then, as for the offline filter, `C^2 + S^2` estimates the (two-sided) power spectral density of the
  signal at `omega_n`, normalised so that the variance is `∫ S(omega) domega / 2π`.
- `:unit_gain`: `A_n = alpha sqrt(alpha^2 + 4omega_n^2) / sqrt(alpha^2 + omega_n^2)`, so that `C` passes a
  signal at exactly `omega_n` with its amplitude preserved, and `S` gives its quadrature with gain
  `omega_n / sqrt(alpha^2 + omega_n^2)`. Then `sqrt(C^2 + S^2)` is its envelope and `atan(S, C)` its
  instantaneous phase (`omega_n*t` plus the phase of the signal). Unlike the offline filter, the one-sided
  window shifts the phases slightly: `C` lags the signal and `S` leads its quadrature, each by about
  `alpha / (2omega_n)` radians, so the envelope and phase have a ripple of relative size about
  `alpha / (2omega_n)` at twice the frequency (for a ripple below 1%, use `alpha ≲ 0.02 omega_n`). At
  `omega_n = 0` this reduces to `alpha`, the single exponential of `set_online_BW_filter_params(N = 1)`.

The `S` terms can be left out of the outputs by setting their coefficients `b_n` to zero.

Arguments
=========
- `alpha`: Decay rate of the exponential window. Must be positive.
- `freqs`: Frequencies (radians per unit time) at which to extract the signal.
- `normalisation`: `:spectral` (default) or `:unit_gain`.

Returns
=======
- A `NamedTuple` of coefficients, `N_coeffs` and `outputs = :separate`.
"""
function set_online_spectrum_filter_params(; alpha::Real, freqs::AbstractVector, normalisation::Symbol = :spectral)
    amplitude(omega) = normalisation === :spectral ? sqrt(2alpha) :
                       alpha * sqrt(alpha^2 + 4omega^2) / sqrt(alpha^2 + omega^2) # Unit gain at omega
    return spectrum_filter_coefficients(amplitude; alpha, freqs, normalisation)
end

# Coefficients of a spectral filter with an exponential window decaying at rate alpha, with a_n = b_n = amplitude(omega_n),
# c_n = alpha and d_n = omega_n at each frequency omega_n in freqs, whose terms are output separately
function spectrum_filter_coefficients(amplitude::Function; alpha::Real, freqs::AbstractVector, normalisation::Symbol)
    alpha > 0 || error("alpha must be positive.")
    length(freqs) > 0 || error("freqs must contain at least one frequency.")
    normalisation in (:spectral, :unit_gain) || error("normalisation must be :spectral or :unit_gain, got $(repr(normalisation))")

    filter_params = NamedTuple()
    for (n, omega) in enumerate(freqs)
        A = amplitude(omega)
        coefficients = NamedTuple{(Symbol("a$n"), Symbol("b$n"), Symbol("c$n"), Symbol("d$n"))}((A, A, alpha, omega))
        filter_params = merge(filter_params, coefficients)
    end
    return merge(filter_params, (; N_coeffs = length(freqs), outputs = :separate))
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
    vars = Dict()
    original_data_filename = config.original_data_filename

    for var_name in var_names_to_filter
        # Only the first frame is needed, so avoid loading the whole time series into memory. The grid and
        # boundary conditions are not read from the file, since those saved by older Oceananigans versions may
        # not be readable (see `_build_templates`); the model uses default boundary conditions.
        fts_data = FieldTimeSeries(original_data_filename, var_name; grid, architecture, backend = InMemory(2), boundary_conditions = nothing)[1]
        vars[Symbol(var_name)] = fts_data
    end
    return NamedTuple(vars)
end

"""
    resolve_map_options(; compute_maps, regrid_to_mean, map_to_mean, lagrangian, normalised, rectilinear, combined_outputs)

Check and reconcile the options controlling the displacement maps ``\\vb*{\\xi}``, for both filter configs:

- `compute_maps`: solve for and output the maps.
- `regrid_to_mean`: interpolate the filtered fields to the mean position, which requires the maps.
- `compute_mean_velocities` (not changed here): output the mean velocities, computed from the maps.

The maps are only displacements from the mean position for a Lagrangian filter (`lagrangian`) with normalised
filter coefficients (`normalised`), so otherwise `regrid_to_mean` is set to `false`. Regridding also currently
requires a rectilinear grid (`rectilinear`), and filtered fields that are output combined rather than as separate
terms (`combined_outputs`, see `filter_output_options`). `map_to_mean` is the name of a removed option, and throws
an error if given.

Returns the reconciled `(compute_maps, regrid_to_mean)`.
"""
function resolve_map_options(; compute_maps, regrid_to_mean, map_to_mean, lagrangian, normalised, rectilinear, combined_outputs)
    if !isnothing(map_to_mean)
        error("The option `map_to_mean` has been removed. Use `compute_maps` to solve for and output the " *
              "displacement maps, and `regrid_to_mean` to interpolate the filtered fields to the mean position.")
    end

    if !lagrangian
        compute_maps && @warn "Advection scheme is `nothing` (Eulerian filter), so the maps are not displacements from the mean position."
        regrid_to_mean && @warn "Advection scheme is `nothing` (Eulerian filter), so setting regrid_to_mean = false."
        regrid_to_mean = false
    end

    if !normalised
        compute_maps && @warn "Filter coefficients are not normalised, so the maps are not displacements from the mean position."
        regrid_to_mean && @warn "Filter coefficients are not normalised, so setting regrid_to_mean = false."
        regrid_to_mean = false
    end

    if !rectilinear && regrid_to_mean
        @warn "The final interpolation to mean position currently only works for RectilinearGrids, so setting regrid_to_mean = false."
        regrid_to_mean = false
    end

    if !combined_outputs && regrid_to_mean
        @warn "The filter outputs its terms separately (outputs = :separate), so there are no combined fields to regrid. Setting regrid_to_mean = false."
        regrid_to_mean = false
    end

    # Regridding uses the maps
    if regrid_to_mean && !compute_maps
        @warn "regrid_to_mean = true requires the maps, so setting compute_maps = true."
        compute_maps = true
    end

    return compute_maps, regrid_to_mean
end

"""
    maps_needed(config::AbstractConfig)

Whether the map variables need to be solved for: they are output if `compute_maps`, and are used to compute
the mean velocities if `compute_mean_velocities`. (`regrid_to_mean` implies `compute_maps`, see `resolve_map_options`.)
"""
maps_needed(config::AbstractConfig) = config.compute_maps || config.compute_mean_velocities

# The number of filter coefficients (a1, b1, c1, d1, a2, ...) in filter_params, not counting its other fields
number_of_coefficient_fields(filter_params::NamedTuple) = count(k -> occursin(r"^[abcd][0-9]+$", String(k)), keys(filter_params))

"""
    create_filtered_vars(config::AbstractConfig)

Creates a `Tuple` of `Symbol`s representing the names of the filtered tracer
variables.

* For a single-exponential filter (`N_coeffs = 0.5`), the function generates
  names with a `_C1` suffix.
* For a multi-coefficient filter (`N_coeffs > 0.5`), it generates pairs of
  names for each coefficient, suffixed with `_C#` and `_S#`, where `#` is the
  coefficient index.

If `compute_maps` or `compute_mean_velocities` is enabled in the configuration,
additional symbols are created for the spatial mapping variables corresponding
to each velocity component, prefixed with `xi_` and suffixed with the corresponding
coefficient names.

Arguments
=========
- `config`: An instance of `AbstractConfig` containing the names of the variables
  to filter, the filter parameters, and the `compute_maps` and `compute_mean_velocities` booleans.

Returns
=======
A `Tuple` of `Symbol`s representing the names of the filtered variables to be
used as tracers in the simulation.
"""
function create_filtered_vars(config::AbstractConfig)
    
    var_names_to_filter = config.var_names_to_filter
    velocity_names = config.velocity_names
    filter_params = config.filter_params
    compute_mean_velocities = config.compute_mean_velocities
    N_coeffs = filter_params.N_coeffs
    label = config.label

    if N_coeffs == 0.5 # special case, single exponential only has a cosine component
        gC_symbols = Symbol[]
        for var_name in var_names_to_filter
            push!(gC_symbols, Symbol(var_name, label, "_C1"))
        end
        # May also need xi maps. We need one for every velocity dimension, so lets use the velocity names to name them
        if maps_needed(config)
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
        if maps_needed(config)
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
    forcing_func = (args...) -> -args[end][1]/(args[end][1]^2 + args[end][2]^2)*args[end-1] # Parameters are last argument, and field dependence is second to last argument. 
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
    _make_gC_relaxation(i::Int, labelled_var_name::String, original_var_name::String, filter_params::NamedTuple, relax_timescale::Real, mask_func::Function, mask_params::Union{NamedTuple, Nothing})
Create a relaxation term for the cosine component (gC) of a filtered variable. The function handles the special case of a single exponential filter where `d` is zero.  

# Arguments
- `i::Int`: The index of the coefficient pair (cᵢ, dᵢ) to use from `filter_params`.
- `labelled_var_name::String`: The name of the variable being filtered (e.g., "T")
    including label if used.
- `original_var_name::String`: The name of the original variable (e.g., "T").
- `filter_params::NamedTuple`: A `NamedTuple` containing all filter coefficients.
- `relax_timescale::Real`: The timescale over which the relaxation occurs.
- `mask_func::Function`: A function that defines the spatial mask for the relaxation.
- `mask_params::Union{NamedTuple, Nothing}`: Optional parameters for the mask function. 

# Returns
- A `Forcing` object configured to compute the relaxation term for the gC field.
"""
function _make_gC_relaxation(i::Int, labelled_var_name::String, original_var_name::String, filter_params::NamedTuple, relax_timescale::Real, mask_func::Function, mask_params::Union{NamedTuple, Nothing})
    if filter_params.N_coeffs == 0.5 # Single exponential special case has a simpler forcing
        c = getproperty(filter_params, Symbol("c",i))
        gCkey = Symbol(labelled_var_name,"_C",i)
        var_key = Symbol(original_var_name)
        # args are (spatial variables, t, field deps, parameters). Parameters are args[end] = (c, relax_timescale, mask_params), and field deps are original variable (args[end-2]), gC (args[end-1])
        # use this general call signature to account for different numbers of spatial variables
        # mask_func takes spatial variables (args[1:end-4] - not time args[end-3]) and mask_params (args[end][3])
        gC_relaxation_func = (args...) -> -1/args[end][2]*(args[end-1] - args[end-2]/args[end][1]) * mask_func(args[1:end-4]...,args[end][3]) 
        return Forcing(gC_relaxation_func, parameters = (c, relax_timescale, mask_params), field_dependencies = (var_key, gCkey))
    else
        c = getproperty(filter_params, Symbol("c",i))
        d = getproperty(filter_params, Symbol("d",i))
        gCkey = Symbol(labelled_var_name, "_C",i)
        var_key = Symbol(original_var_name)
        # args are (spatial variables, t, field deps, parameters). Parameters are args[end] = (c, d, relax_timescale, mask_params), and field deps are original variable (args[end-2]), gC (args[end-1])
        # use this general call signature to account for different numbers of spatial variables
        # mask_func takes spatial variables (args[1:end-4] - not time args[end-3]) and mask_params (args[end][4])
        gC_relaxation_func = (args...) -> -1/args[end][3]*(args[end-1] - args[end-2] * args[end][1]/(args[end][1]^2 + args[end][2]^2)) * mask_func(args[1:end-4]...,args[end][4]) 
        return Forcing(gC_relaxation_func, parameters = (c, d, relax_timescale, mask_params), field_dependencies = (var_key,gCkey))
    end
end

"""
    _make_gS_relaxation(i::Int, labelled_var_name::String, original_var_name::String, filter_params::NamedTuple, relax_timescale::Real, mask_func::Function, mask_params::Union{NamedTuple, Nothing})
Create a relaxation term for the sine component (gS) of a filtered variable.
    
# Arguments
- `i::Int`: The index of the coefficient pair (cᵢ, dᵢ) to use from `filter_params`.
- `labelled_var_name::String`: The name of the variable being filtered (e.g., "T")
    including label if used.
- `original_var_name::String`: The name of the original variable (e.g., "T").
- `filter_params::NamedTuple`: A `NamedTuple` containing all filter coefficients.
- `relax_timescale::Real`: The timescale over which the relaxation occurs.
- `mask_func::Function`: A function that defines the spatial mask for the relaxation.
- `mask_params::Union{NamedTuple, Nothing}`: Optional parameters for the mask function. 

# Returns
- A `Forcing` object configured to compute the relaxation term for the gS field.

"""
function _make_gS_relaxation(i::Int, labelled_var_name::String, original_var_name::String, filter_params::NamedTuple, relax_timescale::Real, mask_func::Function, mask_params::Union{NamedTuple, Nothing})
    c = getproperty(filter_params, Symbol("c",i))
    d = getproperty(filter_params, Symbol("d",i))
    gSkey = Symbol(labelled_var_name, "_S",i)
    var_key = Symbol(original_var_name)
    # args are (spatial variables, t, field deps, parameters). Parameters are args[end] = (c, d, relax_timescale, mask_params), and field deps are original variable (args[end-2]), gS (args[end-1])
    # use this general call signature to account for different numbers of spatial variables
    # mask_func takes spatial variables (args[1:end-4] - not time args[end-3]) and mask_params (args[end][4])
    gS_relaxation_func = (args...) -> -1/args[end][3]*(args[end-1] - args[end-2]*args[end][2]/(args[end][1]^2 + args[end][2]^2)) * mask_func(args[1:end-4]...,args[end][4]) 
    return Forcing(gS_relaxation_func, parameters = (c, d, relax_timescale, mask_params), field_dependencies = (var_key,gSkey))
end

"""
    _make_xiC_relaxation(i::Int, labelled_var_name::String, vel_name::String, filter_params::NamedTuple, relax_timescale::Real, mask_func::Function, mask_params::Union{NamedTuple, Nothing})
Create a relaxation term for the cosine component (xiC) of a map variable. The function handles the special case of a single exponential filter where `d` is zero.  

# Arguments
- `i::Int`: The index of the coefficient pair (cᵢ, dᵢ) to use from `filter_params`.
- `labelled_var_name::String`: The name of the variable being filtered (e.g., "xi_u")
    including label if used.
- `vel_name::String`: The name of the velocity variable (e.g., "u"). 
- `filter_params::NamedTuple`: A `NamedTuple` containing all filter coefficients.
- `relax_timescale::Real`: The timescale over which the relaxation occurs.
- `mask_func::Function`: A function that defines the spatial mask for the relaxation.
- `mask_params::Union{NamedTuple, Nothing}`: Optional parameters for the mask function. 
"""
function _make_xiC_relaxation(i::Int, labelled_var_name::String, vel_name::String, filter_params::NamedTuple, relax_timescale::Real, mask_func::Function, mask_params::Union{NamedTuple, Nothing})
    if filter_params.N_coeffs == 0.5 # Single exponential special case has a simpler forcing
        c = getproperty(filter_params, Symbol("c",i))
        xiCkey = Symbol(labelled_var_name,"_C",i)
        vel_key = Symbol(vel_name)
        # args are (spatial variables, t, field deps, parameters). Parameters are args[end] = (c, relax_timescale, mask_params), and field deps are vel (args[end-2]), xiC (args[end-1])
        # use this general call signature to account for different numbers of spatial variables
        # mask_func takes spatial variables (args[1:end-4] - not time args[end-3]) and mask_params (args[end][3])
        xiC_relaxation_func = (args...) -> -1/args[end][2]*(args[end-1] - (-1/args[end][1]^2)*args[end-2]) * mask_func(args[1:end-4]...,args[end][3]) 
        return Forcing(xiC_relaxation_func, parameters = (c, relax_timescale, mask_params), field_dependencies = (vel_key, xiCkey))
    else
        c = getproperty(filter_params, Symbol("c",i))
        d = getproperty(filter_params, Symbol("d",i))
        xiCkey = Symbol(labelled_var_name, "_C",i)
        vel_key = Symbol(vel_name)
        # args are (spatial variables, t, field deps, parameters). Parameters are args[end] = (c, d, relax_timescale, mask_params), and field deps are velocity (args[end-2]), xiC (args[end-1])
        # use this general call signature to account for different numbers of spatial variables
        # mask_func takes spatial variables (args[1:end-4] - not time args[end-3]) and mask_params (args[end][4])
        xiC_relaxation_func = (args...) -> -1/args[end][3]*(args[end-1] - args[end-2]*(args[end][2]^2 - args[end][1]^2)/(args[end][1]^2 + args[end][2]^2)^2) * mask_func(args[1:end-4]...,args[end][4]) 
        return Forcing(xiC_relaxation_func, parameters = (c, d, relax_timescale, mask_params), field_dependencies = (vel_key, xiCkey))
    end
end


"""
    _make_xiS_relaxation(i::Int, labelled_var_name::String, vel_name::String, filter_params::NamedTuple, relax_timescale::Real, mask_func::Function, mask_params::Union{NamedTuple, Nothing})
Create a relaxation term for the sine component (xiS) of a map variable.

# Arguments
- `i::Int`: The index of the coefficient pair (cᵢ, dᵢ) to use from `filter_params`.
- `labelled_var_name::String`: The name of the variable being filtered (e.g., "xi_u")
    including label if used.
- `vel_name::String`: The name of the velocity variable (e.g., "u"). 
- `filter_params::NamedTuple`: A `NamedTuple` containing all filter coefficients.
- `relax_timescale::Real`: The timescale over which the relaxation occurs.
- `mask_func::Function`: A function that defines the spatial mask for the relaxation.
- `mask_params::Union{NamedTuple, Nothing}`: Optional parameters for the mask function
"""
function _make_xiS_relaxation(i::Int, labelled_var_name::String, vel_name::String, filter_params::NamedTuple, relax_timescale::Real, mask_func::Function, mask_params::Union{NamedTuple, Nothing})
    c = getproperty(filter_params, Symbol("c",i))
    d = getproperty(filter_params, Symbol("d",i))
    xiSkey = Symbol(labelled_var_name, "_S",i)
    vel_key = Symbol(vel_name)
    # args are (spatial variables, t, field deps, parameters). Parameters are args[end] = (c, d, relax_timescale, mask_params), and field deps are velocity (args[end-2]), xiS (args[end-1])
    # use this general call signature to account for different numbers of spatial variables
    # mask_func takes spatial variables (args[1:end-4] - not time args[end-3]) and mask_params (args[end][4])
    xiS_relaxation_func = (args...) -> -1/args[end][3]*(args[end-1] - args[end-2]*(-2*args[end][1]*args[end][2])/(args[end][1]^2 + args[end][2]^2)^2) * mask_func(args[1:end-4]...,args[end][4])
    return Forcing(xiS_relaxation_func, parameters = (c, d, relax_timescale, mask_params), field_dependencies = (vel_key, xiSkey))
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
* For spatial mapping variables (if `compute_maps` or `compute_mean_velocities` is true), the forcing
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
            gC_original_var_forcing = Forcing(original_var_forcing_func, field_dependencies = (;var_key))

            if config.boundary_relaxation
                relax_timescale = config.relax_timescale
                mask_func = config.mask_func
                mask_params = config.mask_params
                gC_relaxation = _make_gC_relaxation(1, labelled_var_name, var_name, filter_params, relax_timescale, mask_func, mask_params)
                gC_forcings_dict[gCkey] = (gC_forcing, gC_original_var_forcing, gC_relaxation)
            else
                gC_forcings_dict[gCkey] = (gC_forcing, gC_original_var_forcing)
            end
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

                if config.boundary_relaxation
                    relax_timescale = config.relax_timescale
                    mask_func = config.mask_func
                    mask_params = config.mask_params
                    xiC_relaxation = _make_xiC_relaxation(1, labelled_var_name, vel_name, filter_params, relax_timescale, mask_func, mask_params)
                    gC_forcings_dict[gCkey] = (xiC_forcing, gC_forcing, xiC_relaxation)
                else
                    gC_forcings_dict[gCkey] = (xiC_forcing, gC_forcing)
                end
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
                gS_forcing_i = _make_gS_forcing(i, labelled_var_name, filter_params)

                if config.boundary_relaxation
                    relax_timescale = config.relax_timescale
                    mask_func = config.mask_func
                    mask_params = config.mask_params
                    gC_relaxation = _make_gC_relaxation(i, labelled_var_name, var_name, filter_params, relax_timescale, mask_func, mask_params)
                    gS_relaxation = _make_gS_relaxation(i, labelled_var_name, var_name, filter_params, relax_timescale, mask_func, mask_params)

                    gC_forcings_dict[gCkey] = (gC_forcing_i, gC_original_var_forcing, gC_relaxation)
                    gS_forcings_dict[gSkey] = (gS_forcing_i, gS_relaxation)
                else
                    gC_forcings_dict[gCkey] = (gC_forcing_i, gC_original_var_forcing)
                    # The forcing for gS is just a filter forcing term
                    gS_forcings_dict[gSkey] = gS_forcing_i
                end
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
                    
                    # The forcing for xiS also includes a term involving xiC and xiS and a term involving the corresponding velocity
                    gS_forcing_i = _make_gS_forcing(i, labelled_var_name, filter_params)
                    xiS_forcing_i = _make_xiS_forcing(i, vel_name, filter_params)

                    if config.boundary_relaxation
                        relax_timescale = config.relax_timescale
                        mask_func = config.mask_func
                        mask_params = config.mask_params
                        xiC_relaxation = _make_xiC_relaxation(i, labelled_var_name, vel_name, filter_params, relax_timescale, mask_func, mask_params)
                        xiS_relaxation = _make_xiS_relaxation(i, labelled_var_name, vel_name, filter_params, relax_timescale, mask_func, mask_params)
                        gC_forcings_dict[gCkey] = (xiC_forcing_i, gC_forcing_i, xiC_relaxation)
                        gS_forcings_dict[gSkey] = (xiS_forcing_i, gS_forcing_i, xiS_relaxation)
                    else
                        gC_forcings_dict[gCkey] = (xiC_forcing_i, gC_forcing_i)
                        gS_forcings_dict[gSkey] = (xiS_forcing_i, gS_forcing_i)
                    end

                end
            end
        end

        return (; NamedTuple(gC_forcings_dict)..., NamedTuple(gS_forcings_dict)...)
    end
end

"""
    create_output_fields(model::AbstractModel, config::AbstractConfig; direction::Symbol = :forward)

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

By default (`outputs = :combined` in `filter_params`), each filtered quantity is output as a
single field. With `outputs = :separate`, each (cosine or sine) term of the filter is output
separately, named after its tracer with a `_scaled` suffix, e.g. `b_C1_scaled`, `xi_u_C1_scaled`,
and `u_C1_scaled` for the mean velocities. With `sine_parity = :odd`, the sine terms are odd in
time, so their backward-pass contributions are negated. See `filter_output_options`.

For the offline filter, the outputs of each pass are defined as that pass's contribution
to the filtered fields, so that the forward and backward outputs are summed (see
[`sum_forward_backward_contributions!`](@ref)). The backward pass is run with the velocities
negated, so its mean velocity outputs are negated (`direction = :backward`).

Arguments
=========
- `model`: An instance of an `AbstractModel` containing the tracer and
  auxiliary fields.
- `config`: An instance of `AbstractConfig` with the names of the variables,
  velocity components, and filter parameters.

Keyword arguments
=================
- `direction`: `:forward` (default, also used for the online filter) or `:backward`, the
  direction of the offline filter pass that the outputs are for.

Returns
=======
A `Dict` where keys are the names of the output fields (e.g.,
`var_name_Lagrangian_filtered`, `xi_vel_name`, `var_name`) and values are the
corresponding reconstructed `Field`s.
"""
function create_output_fields(model::AbstractModel, config::AbstractConfig; direction::Symbol = :forward)
    direction in (:forward, :backward) || error("direction must be :forward or :backward, got :$direction")

    filter_params = config.filter_params
    outputs, sine_parity = filter_output_options(filter_params)

    # Odd sine terms have weight functions that are odd in time, so change sign in the backward pass
    s = (direction === :backward && sine_parity === :odd) ? -1 : 1

    outputs_dict = Dict()
    for q in filtered_output_quantities(config)
        add_filtered_outputs!(outputs_dict, model.tracers, filter_params, q, s, outputs, direction)
    end

    # We can also add the saved vars for comparison if this is an offline filter, otherwise do this manually 
    if (config isa AbstractOfflineConfig) && config.output_original_data
        for var_name in config.var_names_to_filter
            outputs_dict[var_name] = getproperty(model.auxiliary_fields, Symbol(var_name))
        end
        for vel_name in config.velocity_names
            outputs_dict[vel_name] = getproperty(model.velocities, Symbol(vel_name))
        end
    end
 

    return outputs_dict
end

# The backward pass of the offline filter runs with the velocities negated, so its outputs of quantities that
# change sign under time reversal (the mean velocities) are negated to give its contribution to the filtered field.
time_reversed(output, direction) = direction === :backward ? -output : output

"""
    filter_output_options(filter_params::NamedTuple)

Return the output options `(outputs, sine_parity)` of a filter, with their defaults if not given:

- `outputs`: `:combined` (default), to output each filtered quantity as a single field, or `:separate`, to
  output each (cosine or sine) term of the filter separately, named after its tracer with a `_scaled` suffix.
- `sine_parity`: `:even` (default), if the sine terms of the weight function are `sin(d|t|)`, or `:odd`, if they
  are `sin(dt)`. Odd sine terms change sign in the backward pass of the offline filter.
"""
function filter_output_options(filter_params::NamedTuple)
    outputs     = get(filter_params, :outputs, :combined)
    sine_parity = get(filter_params, :sine_parity, :even)
    outputs in (:combined, :separate) || error("filter_params.outputs must be :combined or :separate, got $(repr(outputs))")
    sine_parity in (:even, :odd)      || error("filter_params.sine_parity must be :even or :odd, got $(repr(sine_parity))")
    return outputs, sine_parity
end

# The filtered quantities that are output, one NamedTuple for each with fields
#   quantity:      :tracer (a filtered variable or map) or :velocity (a mean velocity, computed from the maps)
#   combined_name: the name of the combined output
#   separate_stem: the start of the names of the separate outputs (<separate_stem>_C1_scaled, ...)
#   tracer_stem:   the start of the names of the tracers it is computed from (<tracer_stem>_C1, ...)
function filtered_output_quantities(config::AbstractConfig)
    label = config.label

    # When offline filtering, we can turn off advection to get Eulerian filtered fields
    if (config isa AbstractOfflineConfig) && config.advection === nothing
        filter_identifier = "_Eulerian_filtered"
    else
        filter_identifier = "_Lagrangian_filtered"
    end

    row(quantity, combined_name, separate_stem, tracer_stem) = (; quantity, combined_name, separate_stem, tracer_stem)

    quantities = [row(:tracer, var * label * filter_identifier, var * label, var * label) for var in config.var_names_to_filter]
    if config.compute_maps
        append!(quantities, [row(:tracer, "xi_" * vel * label, "xi_" * vel * label, "xi_" * vel * label) for vel in config.velocity_names])
    end
    if config.compute_mean_velocities
        append!(quantities, [row(:velocity, vel * label * filter_identifier, vel * label, "xi_" * vel * label) for vel in config.velocity_names])
    end
    return quantities
end

# The terms (i, has_cosine, has_sine) of a filter that are output separately: terms with a zero coefficient are skipped
filter_terms(filter_params) = filter_params.N_coeffs == 0.5 ? [(1, filter_params.a1 != 0, false)] :
    [(i, getproperty(filter_params, Symbol("a$i")) != 0, getproperty(filter_params, Symbol("b$i")) != 0) for i in 1:filter_params.N_coeffs]

# The coefficients of each term of a filter that is output separately, as a single-term filter, with the suffix of
# its output names (see filtered_output_names)
function separate_term_params(filter_params)
    _, sine_parity = filter_output_options(filter_params)
    terms = Tuple{String, NamedTuple}[]
    for (i, has_cosine, has_sine) in filter_terms(filter_params)
        a, c = getproperty(filter_params, Symbol("a$i")), getproperty(filter_params, Symbol("c$i"))
        if filter_params.N_coeffs == 0.5
            push!(terms, ("_C1_scaled", (; a1 = a, c1 = c, N_coeffs = 0.5)))
        else
            b, d = getproperty(filter_params, Symbol("b$i")), getproperty(filter_params, Symbol("d$i"))
            has_cosine && push!(terms, ("_C$(i)_scaled", (; a1 = a, b1 = zero(b), c1 = c, d1 = d, N_coeffs = 1, sine_parity)))
            has_sine   && push!(terms, ("_S$(i)_scaled", (; a1 = zero(a), b1 = b, c1 = c, d1 = d, N_coeffs = 1, sine_parity)))
        end
    end
    return terms
end

"""
    filtered_output_names(config::AbstractConfig)

Names of the filtered fields output by [`create_output_fields`](@ref), which are combined by
[`sum_forward_backward_contributions!`](@ref).
"""
function filtered_output_names(config::AbstractConfig)
    outputs, _ = filter_output_options(config.filter_params)
    names = String[]
    for q in filtered_output_quantities(config)
        if outputs === :combined
            push!(names, q.combined_name)
        else
            for (i, has_cosine, has_sine) in filter_terms(config.filter_params)
                has_cosine && push!(names, q.separate_stem * "_C$(i)_scaled")
                has_sine   && push!(names, q.separate_stem * "_S$(i)_scaled")
            end
        end
    end
    return Tuple(names)
end

# Coefficients of the C and S tracers in the cosine and sine terms of component i, for a filtered variable or map
# (quantity = :tracer), or a mean velocity computed from the maps (quantity = :velocity, found by integrating the
# filtered velocity by parts). The sine terms are multiplied by s. For s = 1, adding the cosine and sine coefficients
# gives exactly the same floating point operations as the expressions for the combined outputs.
function term_coefficients(filter_params, i, quantity, s)
    a = getproperty(filter_params, Symbol("a$i"))
    c = getproperty(filter_params, Symbol("c$i"))
    if filter_params.N_coeffs == 0.5 # Single exponential: only a cosine term
        return quantity === :tracer ? ((a, 0), nothing) : ((-a * c, 0), nothing)
    end
    b = getproperty(filter_params, Symbol("b$i"))
    d = getproperty(filter_params, Symbol("d$i"))
    if quantity === :tracer
        return (a, 0), (0, s * b)
    else
        return (-a * c, -a * d), (s * b * d, -(s * b) * c)
    end
end

# The linear combination of the C and S tracers with the given coefficients, leaving out zero coefficients
scaled_tracers((cC, cS), gC, gS) = cS == 0 ? cC * gC : cC == 0 ? cS * gS : cC * gC + cS * gS

# Add the outputs for one filtered quantity `q` (see filtered_output_quantities) to `outputs_dict`, either combined
# or separately for each term
function add_filtered_outputs!(outputs_dict, tracers, filter_params, q, s, outputs, direction)
    N_coeffs = filter_params.N_coeffs
    gC(i) = getproperty(tracers, Symbol(q.tracer_stem, "_C", i))
    gS(i) = getproperty(tracers, Symbol(q.tracer_stem, "_S", i))
    sign_for_time_reversal(output) = q.quantity === :velocity ? time_reversed(output, direction) : output

    if outputs === :combined
        total = nothing
        for i in 1:ceil(Int, N_coeffs)
            cosine, sine = term_coefficients(filter_params, i, q.quantity, s)
            term = isnothing(sine) ? cosine[1] * gC(i) :
                                     (cosine[1] + sine[1]) * gC(i) + (cosine[2] + sine[2]) * gS(i)
            total = isnothing(total) ? term : total + term
        end
        outputs_dict[q.combined_name] = sign_for_time_reversal(total)
    else
        for (i, has_cosine, has_sine) in filter_terms(filter_params)
            cosine, sine = term_coefficients(filter_params, i, q.quantity, s)
            gSi = N_coeffs == 0.5 ? nothing : gS(i)
            has_cosine && (outputs_dict[q.separate_stem * "_C$(i)_scaled"] = sign_for_time_reversal(scaled_tracers(cosine, gC(i), gSi)))
            has_sine   && (outputs_dict[q.separate_stem * "_S$(i)_scaled"] = sign_for_time_reversal(scaled_tracers(sine, gC(i), gSi)))
        end
    end
    return outputs_dict
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

Both the filtered variables and the maps are initialised.

Arguments
=========
- `model`: The `AbstractModel` whose tracers are to be initialized.
- `config`: An instance of `AbstractConfig` with the filter parameters.
"""
function initialise_filtered_vars_from_model(model::AbstractModel, config::AbstractConfig)
    filter_params = config.filter_params
    var_names_to_filter = config.var_names_to_filter
    vel_names = config.velocity_names
    label = config.label
    for var_name in var_names_to_filter
        labelled_var_name = var_name * label
        if filter_params.N_coeffs == 0.5 # Special case of single exponential
            filtered_var_C = Symbol(labelled_var_name,"_C1",)
            c1 = filter_params.c1
            # original data can be tracer or auxiliary field
            if Symbol(var_name) in propertynames(model.tracers)
                original_field = getproperty(model.tracers, Symbol(var_name))
            elseif Symbol(var_name) in propertynames(model.auxiliary_fields)
                original_field = getproperty(model.auxiliary_fields, Symbol(var_name))
            else
                error("Variable $var_name not found in model tracers or auxiliary fields.")
            end
            new_field_C = getproperty(model.tracers, filtered_var_C)
            parent(new_field_C) .= 1/c1*parent(original_field)
        else
            for i in 1:filter_params.N_coeffs
                filtered_var_C = Symbol(labelled_var_name,"_C",i)
                filtered_var_S = Symbol(labelled_var_name,"_S",i)
                ci = getproperty(filter_params,Symbol("c$i"))
                di = getproperty(filter_params,Symbol("d$i"))
                if Symbol(var_name) in propertynames(model.tracers)
                    original_field = getproperty(model.tracers,Symbol(var_name))
                elseif Symbol(var_name) in propertynames(model.auxiliary_fields)
                    original_field = getproperty(model.auxiliary_fields,Symbol(var_name))
                else
                    error("Variable $var_name not found in model tracers or auxiliary fields.")
                end
                new_field_C = getproperty(model.tracers, filtered_var_C)
                new_field_S = getproperty(model.tracers, filtered_var_S)
                parent(new_field_C) .= ci/(ci^2 + di^2)*parent(original_field)
                parent(new_field_S) .= di/(ci^2 + di^2)*parent(original_field)
                
            end
        end
    end

    if maps_needed(config)
        for vel_name in vel_names
            if filter_params.N_coeffs == 0.5 # Special case of single exponential
                filtered_map_C = Symbol("xi_", vel_name, label, "_C1")
                c1 = filter_params.c1
                original_vel = getproperty(model.velocities, Symbol(vel_name))
                original_vel_centred = Field(@at (Center, Center, Center) original_vel)
                map_C = getproperty(model.tracers, filtered_map_C)
                parent(map_C) .= (-1/c1^2)*parent(original_vel_centred) 
            else
                for i in 1:filter_params.N_coeffs
                    filtered_map_C = Symbol("xi_", vel_name, label, "_C",i)
                    filtered_map_S = Symbol("xi_", vel_name, label, "_S",i)
                    ci = getproperty(filter_params,Symbol("c$i"))
                    di = getproperty(filter_params,Symbol("d$i"))
                    original_vel = getproperty(model.velocities, Symbol(vel_name))
                    original_vel_centred = Field(@at (Center, Center, Center) original_vel)
                    map_C = getproperty(model.tracers, filtered_map_C)
                    map_S = getproperty(model.tracers, filtered_map_S)
                    parent(map_C) .= ((di^2 - ci^2)/(ci^2 + di^2)^2)*parent(original_vel_centred) 
                    parent(map_S) .= (-2*ci*di/(ci^2 + di^2)^2)*parent(original_vel_centred) 
                end
            end
        end
    end

end

function change_sign_of_map_variables!(model::AbstractModel, config::AbstractConfig)
    # There are no map variables unless we are regridding to the mean position or computing mean velocities
    maps_needed(config) || return nothing

    vel_names = config.velocity_names
    label = config.label
    for vel_name in vel_names
        if config.filter_params.N_coeffs == 0.5 # Special case of single exponential
            filtered_map_C = Symbol("xi_", vel_name, label, "_C1")
            map_C = getproperty(model.tracers, filtered_map_C)
            parent(map_C) .= -parent(map_C)
        else
            for i in 1:config.filter_params.N_coeffs
                filtered_map_C = Symbol("xi_", vel_name, label, "_C",i)
                filtered_map_S = Symbol("xi_", vel_name, label, "_S",i)
                map_C = getproperty(model.tracers, filtered_map_C)
                map_S = getproperty(model.tracers, filtered_map_S)
                parent(map_C) .= -parent(map_C)
                parent(map_S) .= -parent(map_S)
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
    if maps_needed(config)
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

# ──────────────────────────────────────────────────────────────────────────────
# Reading input data
# Input data is read through a BufferedDataReader, directly from the original file.
# ──────────────────────────────────────────────────────────────────────────────

"""
    update_input_data!(model, reader::BufferedDataReader)

Callback that reads directly from the source file via a two-frame
GPU buffer. Calls `advance_buffer!` to slide the window if needed, then
linearly interpolates into the model's velocity and auxiliary fields.
Designed to be used at callsite `UpdateStateCallsite()`, so that the fields are
updated at each substep of multi-stage time steppers.
"""
function update_input_data!(model::AbstractModel, reader::BufferedDataReader)
    t = model.clock.time
    advance_buffer!(reader, t)
    interpolate_to_model!(model, reader, t)
    return nothing
end

"""
    initialise_filtered_vars_from_data(model, reader::BufferedDataReader, config)

Initialises the model's tracer fields, which represent the components of the
filtered variables, to the (scaled) original data at `sim_t = 0` (i.e. the
physical start/end of the filter interval for forward/backward runs). This
improves the "spin-up" of the filter simulation by providing a good starting
point. The data is interpolated from the two buffered frames.

- For a **single-exponential filter** (`N_coeffs = 0.5`), only the `_C1`
  tracer exists and is initialised.
- For a **multi-coefficient filter** (`N_coeffs > 0.5`), both the `_C` and `_S`
  fields for each coefficient are initialised.

If maps are being computed (`compute_maps` or `compute_mean_velocities`), the
maps are initialised from the velocities in the same way.
"""
function initialise_filtered_vars_from_data(model::AbstractModel,
                                            reader::BufferedDataReader,
                                            config::AbstractConfig)
    filter_params = config.filter_params
    label         = config.label

    # Buffer is already primed at t=0 by create_buffered_reader; safe to call again.
    advance_buffer!(reader, 0.0)

    times  = stored_times(reader.source)
    t_phys = reader.direction == :forward ? reader.T_start : reader.T_end
    t_lo   = times[reader.lo_src_idx]
    t_hi   = times[reader.hi_src_idx]
    α = (t_phys - t_lo) / (t_hi - t_lo)
    β = 1 - α

    lo_frames = reader.frames[reader.lo_slot]
    hi_frames = reader.frames[3 - reader.lo_slot]

    # ── Tracer initialisation ─────────────────────────────────────────────────
    for var_name in reader.var_names
        labelled = var_name * label
        lo_f = lo_frames[Symbol(var_name)]
        hi_f = hi_frames[Symbol(var_name)]

        if filter_params.N_coeffs == 0.5
            c1      = filter_params.c1
            field_C = getproperty(model.tracers, Symbol(labelled, "_C1"))
            interior(field_C) .= (1/c1) .* (β .* interior(lo_f) .+ α .* interior(hi_f))
        else
            for i in 1:filter_params.N_coeffs
                ci      = getproperty(filter_params, Symbol("c", i))
                di      = getproperty(filter_params, Symbol("d", i))
                field_C = getproperty(model.tracers, Symbol(labelled, "_C", i))
                field_S = getproperty(model.tracers, Symbol(labelled, "_S", i))
                val = β .* interior(lo_f) .+ α .* interior(hi_f)
                interior(field_C) .= (ci / (ci^2 + di^2)) .* val
                interior(field_S) .= (di / (ci^2 + di^2)) .* val
            end
        end
    end

    # ── Map (xi) initialisation ───────────────────────────────────────────────
    if maps_needed(config)
        for vel_name in reader.vel_names
            lo_v = lo_frames[Symbol(vel_name)]
            hi_v = hi_frames[Symbol(vel_name)]

            # Interpolate each buffered frame to cell-centre, then linearly combine.
            lo_v_c = Field(@at (Center, Center, Center) lo_v)
            hi_v_c = Field(@at (Center, Center, Center) hi_v)
            compute!(lo_v_c)
            compute!(hi_v_c)

            if filter_params.N_coeffs == 0.5
                c1      = filter_params.c1
                field_C = getproperty(model.tracers, Symbol("xi_", vel_name, label, "_C1"))
                interior(field_C) .= (-1/c1^2) .* (β .* interior(lo_v_c) .+ α .* interior(hi_v_c))
            else
                for i in 1:filter_params.N_coeffs
                    ci      = getproperty(filter_params, Symbol("c", i))
                    di      = getproperty(filter_params, Symbol("d", i))
                    field_C = getproperty(model.tracers, Symbol("xi_", vel_name, label, "_C", i))
                    field_S = getproperty(model.tracers, Symbol("xi_", vel_name, label, "_S", i))
                    val = β .* interior(lo_v_c) .+ α .* interior(hi_v_c)
                    interior(field_C) .= ((di^2 - ci^2) / (ci^2 + di^2)^2) .* val
                    interior(field_S) .= (-2*ci*di / (ci^2 + di^2)^2) .* val
                end
            end
        end
    end

    return nothing
end

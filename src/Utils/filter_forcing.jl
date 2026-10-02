# The forcing of the filtered variables: the right-hand sides of the filter equations

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

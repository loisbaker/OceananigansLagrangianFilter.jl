# The forcing of the filtered variables: the right-hand sides of the filter equations
#
# Each filtered tracer (b_C1, b_S1, ..., xi_u_C1, xi_u_S1, ...) has one DiscreteForcing, which
# Oceananigans calls as func(i, j, k, grid, clock, fields, p) at every cell centre. `fields` is
# merge(velocities, tracers, auxiliary_fields), and `p` holds the field names (as `Val`s, so the
# lookup happens at compile time) and the coefficients of term i:
#
#   tracer, cosine:  d/dt gC  = -c gC - d gS + f              [- R/T_relax (gC - f c/(c^2+d^2))]
#   tracer, sine:    d/dt gS  = -c gS + d gC                  [- R/T_relax (gS - f d/(c^2+d^2))]
#   map, cosine:     d/dt xiC = -c/(c^2+d^2) u - c xiC - d xiS [- R/T_relax (xiC - u (d^2-c^2)/(c^2+d^2)^2)]
#   map, sine:       d/dt xiS = -d/(c^2+d^2) u - c xiS + d xiC [- R/T_relax (xiS - u (-2cd)/(c^2+d^2)^2)]
#
# The single exponential (N_coeffs = 0.5) has d = 0 and no sine fields. The bracketed boundary
# relaxation towards the local equilibrium is only included if `relaxation_mask_func` is given;
# R = relaxation_mask_func(x..., relaxation_params) is the relaxation mask, and
# T_relax = relaxation_params.timescale.
#
# The terms are summed in the same order, and the relaxation targets are written in the same
# form, as the previous implementation (a tuple of ContinuousForcings), so the results are
# bit-for-bit unchanged.

# The value of the field called `name` at (i, j, k)
@inline field_value(fields, ::Val{name}, i, j, k) where name = @inbounds getproperty(fields, name)[i, j, k]

# The value of velocity `name` interpolated to (i, j, k), with the operator interp built by velocity_interpolation
@inline velocity_value(fields, ::Val{name}, interp, i, j, k, grid) where name = @inbounds interp(i, j, k, grid, getproperty(fields, name))

# The operator that interpolates a velocity to cell centres, as ContinuousForcing would use for it
velocity_interpolation(vel_name::String) = interpolation_operator(assumed_field_location(Symbol(vel_name)), (Center(), Center(), Center()))

#####
##### Boundary relaxation
#####

"""
    check_relaxation_options(relaxation_mask_func, relaxation_params, grid)

Check the boundary relaxation options of a filter config. If `relaxation_mask_func` is given, it must take one
argument per non-flat dimension of `grid` plus `relaxation_params`, and `relaxation_params` must contain the
relaxation `timescale`.
"""
function check_relaxation_options(relaxation_mask_func, relaxation_params, grid)
    isnothing(relaxation_mask_func) && return nothing

    # We should check that this function only has one method
    length(methods(relaxation_mask_func)) == 1 || @warn "relaxation_mask_func has multiple methods, that could cause issues"

    # We'll also check the number of args is correct
    num_args = first(methods(relaxation_mask_func)).nargs - 2 # First argument is self, last is relaxation_params
    num_non_flat = count(T -> T !== Flat, topology(grid))
    num_args == num_non_flat || error("relaxation_mask_func has the wrong number of arguments: it should take the $num_non_flat non-flat coordinates and relaxation_params")

    relaxation_params isa NamedTuple && haskey(relaxation_params, :timescale) ||
        error("relaxation_params must contain the relaxation timescale, e.g. relaxation_params = (; timescale = 1hour)")
    return nothing
end

# `relaxation` is either `nothing` (no relaxation) or a NamedTuple (mask_func, params), where params contains the timescale.
# (A NamedTuple, rather than a struct, so Oceananigans adapts it to the GPU along with the other parameters.)
boundary_relaxation(config::AbstractConfig) =
    isnothing(config.relaxation_mask_func) ? nothing : (; mask_func = config.relaxation_mask_func, params = config.relaxation_params)

# mask_func takes the non-flat coordinates of the cell centre and relaxation_params
@inline function relaxation_mask(r, i, j, k, grid)
    X = node(i, j, k, grid, Center(), Center(), Center())
    return r.mask_func(X..., r.params)
end

# Add the relaxation of `q` towards `target` to the tendency G
@inline add_relaxation(G, ::Nothing, i, j, k, grid, q, target) = G
@inline add_relaxation(G, r::NamedTuple, i, j, k, grid, q, target) = G + -1 / r.params.timescale * (q - target) * relaxation_mask(r, i, j, k, grid)

#####
##### Tendencies
#####

@inline function tracer_C_tendency(i, j, k, grid, clock, fields, p)
    gC = field_value(fields, p.C, i, j, k)
    gS = field_value(fields, p.S, i, j, k)
    f  = field_value(fields, p.source, i, j, k)
    G  = (-p.c * gC - p.d * gS) + f
    return add_relaxation(G, p.relaxation, i, j, k, grid, gC, f * p.c / (p.c^2 + p.d^2))
end

@inline function tracer_S_tendency(i, j, k, grid, clock, fields, p)
    gC = field_value(fields, p.C, i, j, k)
    gS = field_value(fields, p.S, i, j, k)
    f  = field_value(fields, p.source, i, j, k)
    G  = -p.c * gS + p.d * gC
    return add_relaxation(G, p.relaxation, i, j, k, grid, gS, f * p.d / (p.c^2 + p.d^2))
end

@inline function exponential_tracer_tendency(i, j, k, grid, clock, fields, p)
    gC = field_value(fields, p.C, i, j, k)
    f  = field_value(fields, p.source, i, j, k)
    G  = -p.c * gC + f
    return add_relaxation(G, p.relaxation, i, j, k, grid, gC, f / p.c)
end

@inline function map_C_tendency(i, j, k, grid, clock, fields, p)
    xiC = field_value(fields, p.C, i, j, k)
    xiS = field_value(fields, p.S, i, j, k)
    u  = velocity_value(fields, p.source, p.interp, i, j, k, grid)
    G  = -p.c / (p.c^2 + p.d^2) * u + (-p.c * xiC - p.d * xiS)
    return add_relaxation(G, p.relaxation, i, j, k, grid, xiC, u * (p.d^2 - p.c^2) / (p.c^2 + p.d^2)^2)
end

@inline function map_S_tendency(i, j, k, grid, clock, fields, p)
    xiC = field_value(fields, p.C, i, j, k)
    xiS = field_value(fields, p.S, i, j, k)
    u  = velocity_value(fields, p.source, p.interp, i, j, k, grid)
    G  = -p.d / (p.c^2 + p.d^2) * u + (-p.c * xiS + p.d * xiC)
    return add_relaxation(G, p.relaxation, i, j, k, grid, xiS, u * (-2 * p.c * p.d) / (p.c^2 + p.d^2)^2)
end

@inline function exponential_map_tendency(i, j, k, grid, clock, fields, p)
    xiC = field_value(fields, p.C, i, j, k)
    u  = velocity_value(fields, p.source, p.interp, i, j, k, grid)
    G  = -p.c / (p.c^2 + 0^2) * u + -p.c * xiC
    return add_relaxation(G, p.relaxation, i, j, k, grid, xiC, (-1 / p.c^2) * u)
end

#####
##### Building the forcing
#####

"""
    create_forcing(filtered_vars::Tuple{Vararg{Symbol}}, config::AbstractConfig)

Create the forcing for the filtered tracers named in `filtered_vars` (from [`create_filtered_vars`](@ref)):
the right-hand sides of the filter equations for each filtered variable, and for the map variables if
they are needed, plus the boundary relaxation if `config.relaxation_mask_func` is given.

Returns a `NamedTuple` with one `DiscreteForcing` for each filtered tracer.
"""
function create_forcing(filtered_vars::Tuple{Vararg{Symbol}}, config::AbstractConfig)
    filter_params = config.filter_params
    label = config.label
    exponential = filter_params.N_coeffs == 0.5
    N_terms = exponential ? 1 : filter_params.N_coeffs
    relaxation = boundary_relaxation(config)

    forcing = Dict{Symbol, DiscreteForcing}()

    for i in 1:N_terms
        c = getproperty(filter_params, Symbol("c", i))
        d = exponential ? 0 : getproperty(filter_params, Symbol("d", i))

        for var_name in config.var_names_to_filter
            C = Symbol(var_name, label, "_C", i)
            S = Symbol(var_name, label, "_S", i)
            source = Val(Symbol(var_name))
            if exponential
                forcing[C] = DiscreteForcing(exponential_tracer_tendency, parameters = (; C = Val(C), source, c, relaxation))
            else
                p = (; C = Val(C), S = Val(S), source, c, d, relaxation)
                forcing[C] = DiscreteForcing(tracer_C_tendency, parameters = p)
                forcing[S] = DiscreteForcing(tracer_S_tendency, parameters = p)
            end
        end

        maps_needed(config) || continue

        for vel_name in config.velocity_names
            C = Symbol("xi_", vel_name, label, "_C", i)
            S = Symbol("xi_", vel_name, label, "_S", i)
            source = Val(Symbol(vel_name))
            interp = velocity_interpolation(vel_name)
            if exponential
                forcing[C] = DiscreteForcing(exponential_map_tendency, parameters = (; C = Val(C), source, interp, c, relaxation))
            else
                p = (; C = Val(C), S = Val(S), source, interp, c, d, relaxation)
                forcing[C] = DiscreteForcing(map_C_tendency, parameters = p)
                forcing[S] = DiscreteForcing(map_S_tendency, parameters = p)
            end
        end
    end

    missing_names = setdiff(filtered_vars, keys(forcing))
    isempty(missing_names) || error("No forcing was created for $missing_names: filtered_vars should come from create_filtered_vars(config)")

    return NamedTuple(forcing)
end

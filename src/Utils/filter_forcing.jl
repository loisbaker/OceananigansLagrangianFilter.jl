# The forcing of the filtered variables: the right-hand sides of the filter equations
#
# Each filtered tracer (b_C1, b_S1, ..., xi_u_C1, xi_u_S1, ...) has one DiscreteForcing, which
# Oceananigans calls as func(i, j, k, grid, clock, fields, p) at every cell centre. `fields` is
# merge(velocities, tracers, auxiliary_fields), and `p` holds the field names (as `Val`s, so the
# lookup happens at compile time) and the coefficients of term i:
#
#   tracer, cosine:  d/dt gC  = M (f - c gC - d gS)                [- R/T_relax (gC - f c/(c^2+d^2))]
#   tracer, sine:    d/dt gS  = M (-c gS + d gC)                   [- R/T_relax (gS - f d/(c^2+d^2))]
#   map, cosine:     d/dt xiC = -c/(c^2+d^2) u + M (-c xiC - d xiS) [- R/T_relax (xiC - u (d^2-c^2)/(M (c^2+d^2)^2))]
#   map, sine:       d/dt xiS = -d/(c^2+d^2) u + M (-c xiS + d xiC) [- R/T_relax (xiS - u (-2cd)/(M (c^2+d^2)^2))]
#
# The single exponential (N_coeffs = 0.5) has d = 0 and no sine fields. The bracketed boundary
# relaxation towards the local equilibrium is only included if `relaxation_mask_func` is given;
# R = relaxation_mask_func(x..., relaxation_params) is the relaxation mask, and
# T_relax = relaxation_params.timescale. M = frequency_scale_func(x..., deps..., frequency_scale_params)
# is the frequency scale, or 1 if frequency_scale_func is not given (see docs/src/theory/frequency_scale.md).
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
##### Frequency scale
#####

"""
    check_frequency_scale_options(frequency_scale_func, frequency_scale_params, frequency_scale_dependencies, grid)

Check the frequency scale options of a filter config. If `frequency_scale_func` is given, it must take one argument
per non-flat dimension of `grid`, one per dependency, and `frequency_scale_params`. Without dependencies, it is also
evaluated at every cell centre, and must be finite and positive.
"""
function check_frequency_scale_options(frequency_scale_func, frequency_scale_params, frequency_scale_dependencies, grid)
    if isnothing(frequency_scale_func)
        isempty(frequency_scale_dependencies) || error("frequency_scale_dependencies are given, but there is no frequency_scale_func")
        return nothing
    end

    # We should check that this function only has one method
    length(methods(frequency_scale_func)) == 1 || @warn "frequency_scale_func has multiple methods, that could cause issues"

    # We'll also check the number of args is correct
    num_args = first(methods(frequency_scale_func)).nargs - 2 # First argument is self, last is frequency_scale_params
    num_non_flat = count(T -> T !== Flat, topology(grid))
    num_deps = length(frequency_scale_dependencies)
    num_args == num_non_flat + num_deps ||
        error("frequency_scale_func has the wrong number of arguments: it should take the $num_non_flat non-flat coordinates, the $num_deps frequency_scale_dependencies and frequency_scale_params")

    # With dependencies, M is only known during the run
    if isempty(frequency_scale_dependencies)
        M = CenterField(on_architecture(CPU(), grid))
        set!(M, (X...) -> frequency_scale_func(X..., frequency_scale_params))
        all(m -> isfinite(m) && m > 0, interior(M)) || error("frequency_scale_func must be finite and positive at every cell centre")
    end
    return nothing
end

# `scale` is either `nothing` (a uniform filter) or a NamedTuple (func, dependencies, params). The dependencies are
# `Val`s of the names of fields at cell centres.
function frequency_scale(config::AbstractConfig)
    isnothing(config.frequency_scale_func) && return nothing
    return (; func = config.frequency_scale_func, dependencies = Val.(config.frequency_scale_dependencies),
              params = config.frequency_scale_params)
end

# The frequency scale M at the cell centre (i, j, k), or `nothing` for a uniform filter
@inline frequency_scale_value(::Nothing, i, j, k, grid, fields) = nothing
@inline function frequency_scale_value(s::NamedTuple, i, j, k, grid, fields)
    X = node(i, j, k, grid, Center(), Center(), Center())
    deps = map(name -> field_value(fields, name, i, j, k), s.dependencies)
    return s.func(X..., deps..., s.params)
end

# The same, with the arguments in the order a KernelFunctionOperation passes them
@inline frequency_scale_kernel(i, j, k, grid, scale, fields) = frequency_scale_value(scale, i, j, k, grid, fields)

"""
    check_frequency_scale_dependencies(model, config)

Check that the fields named in `config.frequency_scale_dependencies` are fields of `model` at cell centres, and
not its velocities.
"""
function check_frequency_scale_dependencies(model, config::AbstractConfig)
    fields = merge(model.velocities, model.tracers, model.auxiliary_fields)
    for name in config.frequency_scale_dependencies
        name in keys(fields) || error("frequency_scale_dependencies includes $name, which is not a field of the model")
        name in keys(model.velocities) &&
            error("frequency_scale_dependencies can't include the velocity $name: compute what the frequency scale needs " *
                  "(e.g. the vorticity) as a field at cell centres instead")
        loc = location(getproperty(fields, name))
        loc == (Center, Center, Center) || error("frequency_scale_dependencies must be at cell centres, but $name is at $loc")
    end
    return nothing
end

# The frequency scale M at the cell centres of `model`, as an operation on its current fields, or `nothing` for a
# uniform filter
function frequency_scale_operation(model, config::AbstractConfig)
    scale = frequency_scale(config)
    isnothing(scale) && return nothing
    fields = merge(model.velocities, model.tracers, model.auxiliary_fields)
    return KernelFunctionOperation{Center, Center, Center}(frequency_scale_kernel, model.grid, scale, fields)
end

"""
    frequency_scale_field(model, config)

The frequency scale M at every cell centre of `model`, from its current fields, or `nothing` for a uniform filter.
Checks the dependencies (see [`check_frequency_scale_dependencies`](@ref)), and that M is finite and positive.
"""
function frequency_scale_field(model, config::AbstractConfig)
    isnothing(config.frequency_scale_func) && return nothing
    check_frequency_scale_dependencies(model, config)

    M = Field(frequency_scale_operation(model, config))
    compute!(M)
    all(m -> isfinite(m) && m > 0, Array(interior(M))) ||
        error("frequency_scale_func must be finite and positive, but it is not everywhere for the initial fields")
    return M
end

# Multiply or divide by M, or leave x unchanged for a uniform filter (so that its results are bit-for-bit unchanged)
@inline multiply_by_scale(x, ::Nothing) = x
@inline multiply_by_scale(x, M) = M * x
@inline divide_by_scale(x, ::Nothing) = x
@inline divide_by_scale(x, M) = x / M

#####
##### Tendencies
#####

@inline function tracer_C_tendency(i, j, k, grid, clock, fields, p)
    gC = field_value(fields, p.C, i, j, k)
    gS = field_value(fields, p.S, i, j, k)
    f  = field_value(fields, p.source, i, j, k)
    M  = frequency_scale_value(p.scale, i, j, k, grid, fields)
    G  = multiply_by_scale((-p.c * gC - p.d * gS) + f, M)
    return add_relaxation(G, p.relaxation, i, j, k, grid, gC, f * p.c / (p.c^2 + p.d^2))
end

@inline function tracer_S_tendency(i, j, k, grid, clock, fields, p)
    gC = field_value(fields, p.C, i, j, k)
    gS = field_value(fields, p.S, i, j, k)
    f  = field_value(fields, p.source, i, j, k)
    M  = frequency_scale_value(p.scale, i, j, k, grid, fields)
    G  = multiply_by_scale(-p.c * gS + p.d * gC, M)
    return add_relaxation(G, p.relaxation, i, j, k, grid, gS, f * p.d / (p.c^2 + p.d^2))
end

@inline function exponential_tracer_tendency(i, j, k, grid, clock, fields, p)
    gC = field_value(fields, p.C, i, j, k)
    f  = field_value(fields, p.source, i, j, k)
    M  = frequency_scale_value(p.scale, i, j, k, grid, fields)
    G  = multiply_by_scale(-p.c * gC + f, M)
    return add_relaxation(G, p.relaxation, i, j, k, grid, gC, f / p.c)
end

@inline function map_C_tendency(i, j, k, grid, clock, fields, p)
    xiC = field_value(fields, p.C, i, j, k)
    xiS = field_value(fields, p.S, i, j, k)
    u  = velocity_value(fields, p.source, p.interp, i, j, k, grid)
    M  = frequency_scale_value(p.scale, i, j, k, grid, fields)
    G  = -p.c / (p.c^2 + p.d^2) * u + multiply_by_scale(-p.c * xiC - p.d * xiS, M)
    return add_relaxation(G, p.relaxation, i, j, k, grid, xiC, divide_by_scale(u * (p.d^2 - p.c^2) / (p.c^2 + p.d^2)^2, M))
end

@inline function map_S_tendency(i, j, k, grid, clock, fields, p)
    xiC = field_value(fields, p.C, i, j, k)
    xiS = field_value(fields, p.S, i, j, k)
    u  = velocity_value(fields, p.source, p.interp, i, j, k, grid)
    M  = frequency_scale_value(p.scale, i, j, k, grid, fields)
    G  = -p.d / (p.c^2 + p.d^2) * u + multiply_by_scale(-p.c * xiS + p.d * xiC, M)
    return add_relaxation(G, p.relaxation, i, j, k, grid, xiS, divide_by_scale(u * (-2 * p.c * p.d) / (p.c^2 + p.d^2)^2, M))
end

@inline function exponential_map_tendency(i, j, k, grid, clock, fields, p)
    xiC = field_value(fields, p.C, i, j, k)
    u  = velocity_value(fields, p.source, p.interp, i, j, k, grid)
    M  = frequency_scale_value(p.scale, i, j, k, grid, fields)
    G  = -p.c / (p.c^2 + 0^2) * u + multiply_by_scale(-p.c * xiC, M)
    return add_relaxation(G, p.relaxation, i, j, k, grid, xiC, divide_by_scale((-1 / p.c^2) * u, M))
end

#####
##### Building the forcing
#####

"""
    create_forcing(filtered_vars::Tuple{Vararg{Symbol}}, config::AbstractConfig)

Create the forcing for the filtered tracers named in `filtered_vars` (from [`create_filtered_vars`](@ref)):
the right-hand sides of the filter equations for each filtered variable, and for the map variables if
they are needed, plus the boundary relaxation if `config.relaxation_mask_func` is given. The filter's
frequencies are multiplied by the frequency scale if `config.frequency_scale_func` is given.

Returns a `NamedTuple` with one `DiscreteForcing` for each filtered tracer.
"""
function create_forcing(filtered_vars::Tuple{Vararg{Symbol}}, config::AbstractConfig)
    filter_params = config.filter_params
    label = config.label
    exponential = filter_params.N_coeffs == 0.5
    N_terms = exponential ? 1 : filter_params.N_coeffs
    relaxation = boundary_relaxation(config)
    scale = frequency_scale(config)

    forcing = Dict{Symbol, DiscreteForcing}()

    for i in 1:N_terms
        c = getproperty(filter_params, Symbol("c", i))
        d = exponential ? 0 : getproperty(filter_params, Symbol("d", i))

        for var_name in config.var_names_to_filter
            C = Symbol(var_name, label, "_C", i)
            S = Symbol(var_name, label, "_S", i)
            source = Val(Symbol(var_name))
            if exponential
                forcing[C] = DiscreteForcing(exponential_tracer_tendency, parameters = (; C = Val(C), source, c, relaxation, scale))
            else
                p = (; C = Val(C), S = Val(S), source, c, d, relaxation, scale)
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
                forcing[C] = DiscreteForcing(exponential_map_tendency, parameters = (; C = Val(C), source, interp, c, relaxation, scale))
            else
                p = (; C = Val(C), S = Val(S), source, interp, c, d, relaxation, scale)
                forcing[C] = DiscreteForcing(map_C_tendency, parameters = p)
                forcing[S] = DiscreteForcing(map_S_tendency, parameters = p)
            end
        end
    end

    missing_names = setdiff(filtered_vars, keys(forcing))
    isempty(missing_names) || error("No forcing was created for $missing_names: filtered_vars should come from create_filtered_vars(config)")

    return NamedTuple(forcing)
end

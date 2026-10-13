# The output fields of a Lagrangian filter

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
negated, so its mean velocity outputs are negated (`direction = :backward`). With a frequency scale, the mean
velocities are the velocity of the mean position (see [Varying cutoff frequency](@ref)).

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

    # The frequency scale, or `nothing` for a uniform filter
    M = frequency_scale_operation(model, config)

    outputs_dict = Dict()
    for q in filtered_output_quantities(config)
        add_filtered_outputs!(outputs_dict, model.tracers, filter_params, q, s, outputs, direction, M)
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

# The mean velocities computed from the maps are the velocities in filter time, (u/M)*, so are multiplied by the
# frequency scale M to give the velocity of the mean position (see docs/src/theory/frequency_scale.md)
velocity_of_mean_position(output, ::Nothing) = output
velocity_of_mean_position(output, M) = M * output

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
function add_filtered_outputs!(outputs_dict, tracers, filter_params, q, s, outputs, direction, M)
    N_coeffs = filter_params.N_coeffs
    gC(i) = getproperty(tracers, Symbol(q.tracer_stem, "_C", i))
    gS(i) = getproperty(tracers, Symbol(q.tracer_stem, "_S", i))
    velocity_output(output) = q.quantity === :velocity ? time_reversed(velocity_of_mean_position(output, M), direction) : output

    if outputs === :combined
        total = nothing
        for i in 1:ceil(Int, N_coeffs)
            cosine, sine = term_coefficients(filter_params, i, q.quantity, s)
            term = isnothing(sine) ? cosine[1] * gC(i) :
                                     (cosine[1] + sine[1]) * gC(i) + (cosine[2] + sine[2]) * gS(i)
            total = isnothing(total) ? term : total + term
        end
        outputs_dict[q.combined_name] = velocity_output(total)
    else
        for (i, has_cosine, has_sine) in filter_terms(filter_params)
            cosine, sine = term_coefficients(filter_params, i, q.quantity, s)
            gSi = N_coeffs == 0.5 ? nothing : gS(i)
            has_cosine && (outputs_dict[q.separate_stem * "_C$(i)_scaled"] = velocity_output(scaled_tracers(cosine, gC(i), gSi)))
            has_sine   && (outputs_dict[q.separate_stem * "_S$(i)_scaled"] = velocity_output(scaled_tracers(sine, gC(i), gSi)))
        end
    end
    return outputs_dict
end

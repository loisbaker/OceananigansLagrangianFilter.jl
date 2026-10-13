# Initialising the filtered variables, and the closure for them

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

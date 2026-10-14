# The original and filtered variables of a Lagrangian filter

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
A `NamedTuple` where each key is a `Symbol` of a variable name to be filtered, or of an
extra auxiliary variable (`config.extra_auxiliary_var_names`), and each value is a field
for that variable.
"""
function create_original_vars(config::AbstractConfig)

    var_names = (config.var_names_to_filter..., config.extra_auxiliary_var_names...)
    grid = config.grid
    architecture = config.architecture
    vars = Dict()
    original_data_filename = config.original_data_filename

    for var_name in var_names
        # Only the first frame is needed, so avoid loading the whole time series into memory. The grid and
        # boundary conditions are not read from the file, since those saved by older Oceananigans versions may
        # not be readable (see `_build_templates`); the model uses default boundary conditions.
        fts_data = FieldTimeSeries(original_data_filename, var_name; grid, architecture, backend = InMemory(2), boundary_conditions = nothing)[1]
        vars[Symbol(var_name)] = fts_data
    end
    return NamedTuple(vars)
end

"""
    resolve_map_options(; compute_maps, regrid_to_mean, map_to_mean, lagrangian, normalised, supported_grid, combined_outputs)

Check and reconcile the options controlling the displacement maps ``\\vb*{\\xi}``, for both filter configs:

- `compute_maps`: solve for and output the maps.
- `regrid_to_mean`: interpolate the filtered fields to the mean position, which requires the maps.
- `compute_mean_velocities` (not changed here): output the mean velocities, computed from the maps.

The maps are only displacements from the mean position for a Lagrangian filter (`lagrangian`) with normalised
filter coefficients (`normalised`), so otherwise `regrid_to_mean` is set to `false`. Regridding also currently
requires a `RectilinearGrid` or `LatitudeLongitudeGrid` (`supported_grid`), and filtered fields that are output combined rather than as separate
terms (`combined_outputs`, see `filter_output_options`). `map_to_mean` is the name of a removed option, and throws
an error if given.

Returns the reconciled `(compute_maps, regrid_to_mean)`.
"""
function resolve_map_options(; compute_maps, regrid_to_mean, map_to_mean, lagrangian, normalised, supported_grid, combined_outputs)
    if !isnothing(map_to_mean)
        error("The option `map_to_mean` has been removed. Use `compute_maps` to solve for and output the " *
              "displacement maps, and `regrid_to_mean` to interpolate the filtered fields to the mean position.")
    end

    if !lagrangian
        compute_maps && @warn "lagrangian = false (Eulerian filter), so the maps are not displacements from the mean position."
        regrid_to_mean && @warn "lagrangian = false (Eulerian filter), so setting regrid_to_mean = false."
        regrid_to_mean = false
    end

    if !normalised
        compute_maps && @warn "Filter coefficients are not normalised, so the maps are not displacements from the mean position."
        regrid_to_mean && @warn "Filter coefficients are not normalised, so setting regrid_to_mean = false."
        regrid_to_mean = false
    end

    if !supported_grid && regrid_to_mean
        @warn "The final interpolation to mean position currently only works for RectilinearGrids and LatitudeLongitudeGrids, so setting regrid_to_mean = false."
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

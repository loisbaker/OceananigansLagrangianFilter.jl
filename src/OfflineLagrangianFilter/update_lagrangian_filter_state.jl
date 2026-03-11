using Oceananigans: UpdateStateCallsite
using Oceananigans.BoundaryConditions
using Oceananigans.BoundaryConditions: update_boundary_conditions! 
using Oceananigans.TurbulenceClosures: compute_diffusivities!
using Oceananigans.Fields: compute!
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using Oceananigans.Models: update_model_field_time_series!

import Oceananigans.TimeSteppers: update_state!

"""
    update_state!(model::LagrangianFilter, callbacks=[])

Update peripheral aspects of the model (halo regions, diffusivities) to the current model state. If `callbacks` are provided (in an array),
they are called in the end.
"""
function update_state!(model::LagrangianFilter, callbacks=[]; compute_tendencies = true)
    
    # Mask immersed tracers
    foreach(model.tracers) do tracer
        mask_immersed_field!(tracer)
    end

    # Update all FieldTimeSeries used in the model
    update_model_field_time_series!(model, model.clock)

    # Update the boundary conditions
    update_boundary_conditions!(fields(model), model)

    # Fill halos for velocities and tracers
    fill_halo_regions!(model.tracers, model.clock, fields(model); fill_open_bcs = false, async = true)

    # Compute auxiliary fields
    for aux_field in model.auxiliary_fields
        compute!(aux_field)
    end

    # Calculate diffusivities 
    compute_auxiliaries!(model)
    fill_halo_regions!(model.closure_fields; only_local_halos = true)
    
    for callback in callbacks
        callback.callsite isa UpdateStateCallsite && callback(model)
    end

    # Keep functionality here to not compute tendencies when set is called if needed. 
    # This allows the `update_input_data!` callback to use `set!`` (which calls `update_state!` without the callback) 
    # without recomputing tendencies each time (for them just to be recomputed here again).
    compute_tendencies && compute_tendencies!(model, callbacks)

    return nothing
end

function compute_auxiliaries!(model::LagrangianFilter; κ_parameters = :xyz) 

    closure = model.closure
    diffusivity = model.closure_fields

    # Compute diffusivities
    compute_diffusivities!(diffusivity, closure, model; parameters = κ_parameters)

    return nothing
end

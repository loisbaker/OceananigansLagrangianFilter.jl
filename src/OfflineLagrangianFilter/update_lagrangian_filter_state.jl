using Oceananigans: UpdateStateCallsite
using Oceananigans.Advection: update_advection_timestep!
using Oceananigans.BoundaryConditions
using Oceananigans.BoundaryConditions: update_boundary_conditions! 
using Oceananigans.TurbulenceClosures: compute_closure_fields!
import Oceananigans.TurbulenceClosures: step_closure_prognostics!
using Oceananigans.Fields: compute!
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using Oceananigans.Models: update_model_field_time_series!

import Oceananigans.TimeSteppers: update_state!

"""
    update_state!(model::LagrangianFilter, callbacks=[])

Update peripheral aspects of the model (halo regions, closure_fields) to the current model state. If `callbacks` are provided (in an array),
they are called in the end.
"""
function update_state!(model::LagrangianFilter, callbacks=[])
    
    # Mask immersed tracers
    foreach(model.tracers) do tracer
        mask_immersed_field!(tracer)
    end

    # Update all FieldTimeSeries used in the model
    update_model_field_time_series!(model, model.clock)

    # Update the boundary conditions (this doesn't do anything by default for most BCs)
    # update_boundary_conditions!(fields(model), model)

    # Fill halos for tracers - we're only applying the BCs to the tracers - not the velocities or the auxiliary fields. 
    fill_halo_regions!(model.tracers, model.clock, fields(model); 
                       fill_normal_flow_bcs = false, async = true)

    # Compute auxiliary fields
    for aux_field in model.auxiliary_fields
        compute!(aux_field)
    end

    # Calculate closure fields 
    compute_auxiliaries!(model)
    fill_halo_regions!(model.closure_fields; only_local_halos = true)
    
    for callback in callbacks
        callback.callsite isa UpdateStateCallsite && callback(model)
    end

    update_advection_timestep!(model.advection, model.timestepper, model.clock)
    compute_tendencies!(model, callbacks)

    return nothing
end


function compute_auxiliaries!(model::LagrangianFilter; κ_parameters = :xyz) 

    closure = model.closure
    closure_fields = model.closure_fields

    # Compute closure fields
    compute_closure_fields!(closure_fields, closure, model; parameters = κ_parameters)

    return nothing
end

step_closure_prognostics!(model::LagrangianFilter, Δt) =
    step_closure_prognostics!(model.closure_fields, model.closure, model, Δt)

using Oceananigans: fields, TendencyCallsite
using Oceananigans.Utils: work_layout
using Oceananigans.Models: complete_communication_and_compute_buffer!, interior_tendency_kernel_parameters
using Oceananigans.Grids: get_active_cells_map

import Oceananigans.TimeSteppers: compute_tendencies!
import Oceananigans.TimeSteppers: compute_flux_bc_tendencies!

"""
    compute_tendencies!(model::NonhydrostaticModel, callbacks)

Calculate the interior and boundary contributions to tendency terms without the
contribution from non-hydrostatic pressure.
"""
function compute_tendencies!(model::LagrangianFilter, callbacks)

    # Note:
    #
    # "tendencies" is a NamedTuple of OffsetArrays corresponding to the tendency data for use
    # in GPU computations.
    #
    # "model.timestepper.Gⁿ" is a NamedTuple of Fields, whose data also corresponds to
    # tendency data.

    grid = model.grid
    arch = architecture(grid)

    # Calculate contributions to momentum and tracer tendencies from fluxes and volume terms in the
    # interior of the domain
    kernel_parameters = interior_tendency_kernel_parameters(arch, grid)
    active_cells_map  = get_active_cells_map(model.grid, Val(:interior))
    
    compute_interior_tendency_contributions!(model, kernel_parameters; active_cells_map)
    complete_communication_and_compute_buffer!(model, grid, arch)

    # # Calculate contributions to momentum and tracer tendencies from user-prescribed fluxes across the
    # # boundaries of the domain
    # compute_boundary_tendency_contributions!(model.timestepper.Gⁿ,
    #                                          model.architecture,
    #                                          model.velocities,
    #                                          model.tracers,
    #                                          model.clock,
    #                                          fields(model))

    for callback in callbacks
        callback.callsite isa TendencyCallsite && callback(model)
    end

    return nothing
end

""" Store previous value of the source term and compute current source term. """
function compute_interior_tendency_contributions!(model, kernel_parameters; active_cells_map = nothing)

    tendencies           = model.timestepper.Gⁿ
    arch                 = model.architecture
    grid                 = model.grid
    advection            = model.advection
    closure              = model.closure
    velocities           = model.velocities
    tracers              = model.tracers
    auxiliary_fields     = model.auxiliary_fields
    diffusivities        = model.closure_fields
    forcings             = model.forcing
    clock                = model.clock

    
    start_tracer_kernel_args = (advection, closure)
    end_tracer_kernel_args   = (velocities, tracers, auxiliary_fields, diffusivities)

    for tracer_index in 1:length(tracers)
        @inbounds c_tendency = tendencies[tracer_index + 3] # This assumes that tracers are after velocities in the NamedTuple, and comes from prognostic_fields
        @inbounds forcing = forcings[tracer_index + 3] # This assumes that tracers are after velocities in the NamedTuple, and comes from prognostic_fields
        @inbounds c_immersed_bc = tracers[tracer_index].boundary_conditions.immersed
        @inbounds tracer_name = keys(tracers)[tracer_index]

        args = tuple(Val(tracer_index), Val(tracer_name),
                     start_tracer_kernel_args..., 
                     c_immersed_bc,
                     end_tracer_kernel_args...,
                     clock, forcing)

        launch!(arch, grid, kernel_parameters, compute_Gc!, 
                c_tendency, grid, args;
                active_cells_map)
    end

    return nothing
end

#####
##### Tracer(s)
#####

""" Calculate the right-hand-side of the tracer advection-diffusion equation. """
@kernel function compute_Gc!(Gc, grid, args)
    i, j, k = @index(Global, NTuple)
    @inbounds Gc[i, j, k] = tracer_tendency(i, j, k, grid, args...)
end


#####
##### Boundary contributions to tendencies due to user-prescribed fluxes
#####

""" Apply boundary conditions by adding flux divergences to the right-hand-side. """
function compute_flux_bc_tendencies!(model::LagrangianFilter)
    
    Gⁿ    = model.timestepper.Gⁿ
    arch  = model.architecture
    clock = model.clock
    model_fields = fields(model)
    prognostic_fields = merge(model.velocities, model.tracers)

    foreach(i -> compute_x_bcs!(Gⁿ[i], prognostic_fields[i], arch, clock, model_fields), 1:length(prognostic_fields))
    foreach(i -> compute_y_bcs!(Gⁿ[i], prognostic_fields[i], arch, clock, model_fields), 1:length(prognostic_fields))
    foreach(i -> compute_z_bcs!(Gⁿ[i], prognostic_fields[i], arch, clock, model_fields), 1:length(prognostic_fields))
    return nothing
end


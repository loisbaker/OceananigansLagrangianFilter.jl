using Oceananigans.TimeSteppers: _ab2_step_field!, implicit_step!
import Oceananigans.TimeSteppers: ab2_step!

"""
    ab2_step!(model::LagrangianFilter, Δt, callbacks)

Advance `LagrangianFilter` by one Adams-Bashforth 2nd-order time step without pressure correction.
Dispatches to `pressure_correction_ab2_step!` which implements a predictor-corrector scheme
"""
ab2_step!(model::LagrangianFilter, args...) =
    pressure_correction_ab2_step!(model, args...)

"""
    pressure_correction_ab2_step!(model, Δt, callbacks)

Implement the AB2 time step with pressure correction for `LagrangianFilter`.

This predictor-corrector scheme is based on the NonhydrostaticModel AB2 step, but
does not advance velocities or correct pressure. The step proceeds as follows:

1. Advances tracers: `cⁿ⁺¹ = cⁿ + Δt * AB2(Gᶜ)`
2. Applies implicit vertical diffusion (if configured)
"""
function pressure_correction_ab2_step!(model, Δt, callbacks)
    grid = model.grid

    # Compute flux bc tendencies
    compute_flux_bc_tendencies!(model)

    # Velocity steps
    # for (i, field) in enumerate(model.velocities)
    #     kernel_args = (field, Δt, model.timestepper.χ, model.timestepper.Gⁿ[i], model.timestepper.G⁻[i])
    #     launch!(architecture(grid), grid, :xyz, _ab2_step_field!, kernel_args...; exclude_periphery=true)

    #     implicit_step!(field,
    #                    model.timestepper.implicit_solver,
    #                    model.closure,
    #                    model.closure_fields,
    #                    nothing,
    #                    model.clock,
    #                    fields(model),
    #                    Δt)
    # end

    # Tracer steps
    for (i, name) in enumerate(propertynames(model.tracers))
        field = model.tracers[name]
        kernel_args = (field, Δt, model.timestepper.χ, model.timestepper.Gⁿ[name], model.timestepper.G⁻[name])
        launch!(architecture(grid), grid, :xyz, _ab2_step_field!, kernel_args...)

        implicit_step!(field,
                       model.timestepper.implicit_solver,
                       model.closure,
                       model.closure_fields,
                       Val(i),
                       model.clock,
                       fields(model),
                       Δt)
    end

    # compute_pressure_correction!(model, Δt)
    # make_pressure_correction!(model, Δt)

    return nothing
end
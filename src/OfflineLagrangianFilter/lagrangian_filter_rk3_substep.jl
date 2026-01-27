using Oceananigans.TimeSteppers: _rk3_substep_field!, stage_Δt
import Oceananigans.TimeSteppers: rk3_substep!

"""
    rk3_substep!(model::LagrangianFilter, Δt, γ, ζ, callbacks)

Perform a single RK3 substep for `LagrangianFilter` tracers without pressure correction.
Dispatches to `tracer_rk3_substep!` which advances tracers
using the RK3 coefficients as in NonhydrostaticModel, but doesn't advance velocities or correct pressure. 
"""
rk3_substep!(model::LagrangianFilter, Δt, γ, ζ, callbacks) =
    tracer_rk3_substep!(model, Δt, γ, ζ, callbacks)

"""
    tracer_rk3_substep!(model, Δt, γⁿ, ζⁿ, callbacks)

Implement a single RK3 substep for `LagrangianFilter`.

The substep advances the state as

    U += Δt * (γⁿ * Gⁿ + ζⁿ * G⁻)

where:
- `γⁿ` is the coefficient for the current tendency
- `ζⁿ` is the coefficient for the previous tendency (or `nothing` for the first substep)
- The effective substep size is `Δτ = Δt * (γⁿ + ζⁿ)`

After advancing velocities, a pressure Poisson equation is solved and velocities
are corrected to satisfy the incompressibility constraint.
"""
function tracer_rk3_substep!(model, Δt, γⁿ, ζⁿ, callbacks)
    Δτ   = stage_Δt(Δt, γⁿ, ζⁿ)
    grid = model.grid

    compute_flux_bc_tendencies!(model)

    # Tracer steps
    for (i, name) in enumerate(propertynames(model.tracers))
        field = model.tracers[name]
        kernel_args = (field, Δt, γⁿ, ζⁿ, model.timestepper.Gⁿ[name], model.timestepper.G⁻[name])
        launch!(architecture(grid), grid, :xyz, _rk3_substep_field!, kernel_args...)

        implicit_step!(field,
                       model.timestepper.implicit_solver,
                       model.closure,
                       model.closure_fields,
                       Val(i),
                       model.clock,
                       fields(model),
                       Δτ)
    end


    return nothing
end
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.TimeSteppers: update_state!

import Oceananigans.Fields: set!

"""
    set!(model::LagrangianFilter; kwargs...)

Set velocity and tracer fields of `model`. The keyword arguments
`kwargs...` take the form `name=data`, where `name` refers to one of the
fields of `model.velocities` or `model.tracers`, and the `data` may be an array,
a function with arguments `(x, y, z)`, or any data type for which a
`set!(ϕ::AbstractField, data)` function exists.

"""
function set!(model::LagrangianFilter; kwargs...)
    for (fldname, value) in kwargs
        if fldname ∈ propertynames(model.velocities)
            ϕ = getproperty(model.velocities, fldname)
        elseif fldname ∈ propertynames(model.tracers)
            ϕ = getproperty(model.tracers, fldname)
        else
            throw(ArgumentError("name $fldname not found in model.velocities or model.tracers."))
        end
        set!(ϕ, value)

        fill_halo_regions!(ϕ, model.clock, fields(model))
    end

    # Apply a mask
    foreach(mask_immersed_field!, model.tracers)
    foreach(mask_immersed_field!, model.velocities)
    update_state!(model; compute_tendencies = false)

    return nothing
end

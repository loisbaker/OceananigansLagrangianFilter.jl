using Oceananigans.Advection
using Oceananigans.Operators

using Oceananigans.TurbulenceClosures: ∂ⱼ_τ₁ⱼ, ∂ⱼ_τ₂ⱼ, ∂ⱼ_τ₃ⱼ, ∇_dot_qᶜ
using Oceananigans.TurbulenceClosures: immersed_∂ⱼ_τ₁ⱼ, immersed_∂ⱼ_τ₂ⱼ, immersed_∂ⱼ_τ₃ⱼ, immersed_∇_dot_qᶜ
using Oceananigans.Forcings: with_advective_forcing



"""
    $(SIGNATURES)

Return the tendency for a tracer field with index `tracer_index`
at grid point `i, j, k`.

The tendency is called ``G_c`` and defined via

```math
∂_t c = G_c ,
```

where `c = C[tracer_index]`.

`closure` and `buoyancy` are types encoding information about the prescribed
turbulence closure and buoyancy model.

`background_fields` is a `NamedTuple` containing background velocity and tracer
`FunctionFields`.

The arguments `velocities`, `tracers`, and `diffusivities` are `NamedTuple`s with the three
velocity components, tracer fields, and precalculated diffusivities where applicable.
`forcings` is a named tuple of forcing functions.

`clock` keeps track of `clock.time` and `clock.iteration`.
"""
@inline function tracer_tendency(i, j, k, grid,
                                 val_tracer_index::Val{tracer_index},
                                 val_tracer_name,
                                 advection,
                                 closure,
                                 c_immersed_bc,
                                 velocities,
                                 tracers,
                                 auxiliary_fields,
                                 diffusivities,
                                 clock,
                                 forcing) where tracer_index

    @inbounds c = tracers[tracer_index]

    model_fields = merge(velocities, tracers, auxiliary_fields)


    total_velocities = with_advective_forcing(forcing, velocities)
    buoyancy = nothing
    return ( - div_Uc(i, j, k, grid, advection, total_velocities, c)
             + c_div_U(i, j, k, grid, advection, total_velocities, c)
             - ∇_dot_qᶜ(i, j, k, grid, closure, diffusivities, val_tracer_index, c, clock, model_fields, buoyancy) 
             - immersed_∇_dot_qᶜ(i, j, k, grid, c, c_immersed_bc, closure, diffusivities, val_tracer_index, clock, model_fields)
             + forcing(i, j, k, grid, clock, model_fields))
end

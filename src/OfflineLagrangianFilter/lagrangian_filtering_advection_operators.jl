using Oceananigans.Fields: ZeroField

#####
##### Tracer advection operator for when velocity is divergent
#####

"""
    c_div_U(i, j, k, grid, advection, U, c)

Calculate the product of the tracer concentration ``c`` with 
the horizontal divergence of the velocity field ``𝐔 = (u, v)``, ``c ∇·𝐔``,

```
c[i, j, k] * divᶜᶜᶜ(U)
```

which ends up at the location `ccc`.
"""
@inline function c_div_U(i, j, k, grid, advection, U, c)
    return c[i, j, k] * divᶜᶜᶜ(i, j, k, grid, U.u, U.v, U.w)
end

const ZeroU = NamedTuple{(:u, :v, :w), Tuple{ZeroField, ZeroField, ZeroField}}

# Fallbacks for zero velocities, zero tracer and `nothing` advection
@inline c_div_U(i, j, k, grid, advection, ::ZeroU, c) = zero(grid)
@inline c_div_U(i, j, k, grid, advection, U, ::ZeroField) = zero(grid)
@inline c_div_U(i, j, k, grid, advection, ::ZeroU, ::ZeroField) = zero(grid)

@inline c_div_U(i, j, k, grid, ::Nothing, U, c) = zero(grid)
@inline c_div_U(i, j, k, grid, ::Nothing, ::ZeroU, c) = zero(grid)
@inline c_div_U(i, j, k, grid, ::Nothing, U, ::ZeroField) = zero(grid)
@inline c_div_U(i, j, k, grid, ::Nothing, ::ZeroU, ::ZeroField) = zero(grid)

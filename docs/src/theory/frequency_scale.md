# Spatially varying cutoff equations

For a positive, stationary mask ``M(\vb*{x})``, define a filter clock along
particle trajectories ``\vb*{\varphi}(\vb*{a},t)`` by

```math
\tau(t)-\tau(s)=\int_s^t M(\vb*{\varphi}(\vb*{a},r))\,\mathrm{d}r.
```

The offline mean uses the reference kernel ``G`` in filter time. Converting the
integral to physical time ``s`` introduces the Jacobian
``\mathrm{d}\tau(s)/\mathrm{d}s=M(\vb*{\varphi}(\vb*{a},s))``:

```math
f^*(\vb*{\varphi}(\vb*{a},t),t)=\int_{-\infty}^{\infty}
G(\tau(t)-\tau(s)) f(\vb*{\varphi}(\vb*{a},s),s)
M(\vb*{\varphi}(\vb*{a},s))\,\mathrm{d}s.
```

The Jacobian makes the kernel integrate to one, preserving constants in the
infinite-window limit. If ``M=m`` along a trajectory, the kernel is
``mG(m(t-s))`` and the local cutoff is ``m\,\mathrm{freq_c}``. If ``M`` varies,
there is no single physical-time frequency response.

With ``D_t=\partial_t+\boldsymbol{u}\cdot\nabla``, each forward filter pair obeys

```math
D_t g_{Ck}=M(f-c_k g_{Ck}-d_k g_{Sk}),\qquad
D_t g_{Sk}=M(-c_k g_{Sk}+d_k g_{Ck}).
```

Reconstruction uses the reference coefficients ``a_k,b_k``. The position-map
equations retain the reference velocity source and multiply their decay terms
by ``M``. Their local equilibrium, also used for initialization and boundary
relaxation, is

```math
\boldsymbol{\xi}_{Ck}
=\frac{d_k^2-c_k^2}{M(c_k^2+d_k^2)^2}\boldsymbol{u},\qquad
\boldsymbol{\xi}_{Sk}
=\frac{-2c_kd_k}{M(c_k^2+d_k^2)^2}\boldsymbol{u}.
```

These equilibria do not account for a particle's prehistory across mask
gradients; finite runs retain endpoint transients. For the single-exponential
case, set ``d_1=0`` and omit the sine variables.

## Online filtering

For a stationary, positive spatial cutoff mask, the online filter uses the
forward part of the [spatially varying cutoff equations](@ref "Spatially varying cutoff equations"):
the filter clock advances at ``\mathrm{d}\tau/\mathrm{d}t=M``, and both the
tracer input and filter decay are multiplied by ``M``.

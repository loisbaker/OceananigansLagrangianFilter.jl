# Varying cutoff frequency

For a positive frequency scale ``M(\vb*{x}, t)`` (`frequency_scale_func`), which multiplies the filter's
frequencies locally, define a filter clock along particle trajectories ``\vb*{\varphi}(\vb*{a},t)`` by

```math
\tau(t)-\tau(s)=\int_s^t M(\vb*{\varphi}(\vb*{a},r), r)\,\mathrm{d}r.
```

The offline mean uses the reference kernel ``G`` in filter time. Converting the
integral to physical time ``s`` introduces the Jacobian
``\mathrm{d}\tau(s)/\mathrm{d}s=M(\vb*{\varphi}(\vb*{a},s), s)``:

```math
f^*(\vb*{\varphi}(\vb*{a},t),t)=\int_{-\infty}^{\infty}
G(\tau(t)-\tau(s)) f(\vb*{\varphi}(\vb*{a},s),s)
M(\vb*{\varphi}(\vb*{a},s), s)\,\mathrm{d}s.
```

The Jacobian makes the kernel integrate to one, preserving constants in the
infinite-window limit. If ``M=m`` along a trajectory, the kernel is
``mG(m(t-s))`` and the local cutoff is ``m\,\mathrm{freq_c}``. If ``M`` varies,
there is no single physical-time frequency response.

With ``D_t=\partial_t+\vb*{u}\cdot\nabla``, each forward filter pair obeys

```math
D_t g_{Ck}=M(f-c_k g_{Ck}-d_k g_{Sk}),\qquad
D_t g_{Sk}=M(-c_k g_{Sk}+d_k g_{Ck}).
```

Reconstruction uses the reference coefficients ``a_k,b_k``. The position-map
equations retain the reference velocity source and multiply their decay terms
by ``M``. Their local equilibrium, also used for initialization and boundary
relaxation, is

```math
\vb*{\xi}_{Ck}
=\frac{d_k^2-c_k^2}{M(c_k^2+d_k^2)^2}\vb*{u},\qquad
\vb*{\xi}_{Sk}
=\frac{-2c_kd_k}{M(c_k^2+d_k^2)^2}\vb*{u}.
```

These equilibria do not account for a particle's prehistory across gradients of
``M``; finite runs retain endpoint transients. For the single-exponential
case, set ``d_1=0`` and omit the sine variables.

## Frequency scales that vary in time

``M`` can also vary in time, through `frequency_scale_dependencies`: the values of other fields (for example
the vorticity) at the same point. The filter clock above is unchanged, with ``M`` evaluated at the particle's
position at each time. Offline, the forward and backward passes see the same ``M(\vb*{x}, t)``, since the fields it
depends on are not reversed in sign on the backward pass, so the two passes use the same filter clock. The local
equilibria used for initialisation and boundary relaxation use the current value of ``M``.

``M`` must stay positive: it is the rate at which the filter clock runs, so if it reached zero the filter would
stop, and if it became negative the clock would run backwards and the filter would amplify rather than damp. With
dependencies, ``M`` is only checked when the filtered variables are initialised, so choose `frequency_scale_func`
so that it is positive for any values the fields can take, for example ``M = 1 + |\zeta|/f`` rather than
``M = 1 + \zeta/f``.

Offline, ``M`` can't depend on the filtered fields themselves. The offline filter is the sum of a forward pass,
for the past, and a backward pass, for the future, and both must use the same ``M`` along each trajectory for the
two halves to make up one filter with one filter clock. During each pass the filter's tracers only hold that pass's
contribution, so an ``M`` computed from them would differ between the passes. Instead, filter the field first
(e.g. with a uniform cutoff), save it to the data file, and read it in with `extra_auxiliary_var_names` for a
second filter run. Online, the filtered fields exist during the run, so ``M`` can depend on them, through an
auxiliary field kept equal to the filtered field; ``M`` then depends on the filter's own output, so the filter is
no longer linear.

Whether ``M`` varies in space or time, the filter has no single frequency response or time shift, so
[`compute_time_shift!`](@ref) is not available. The Eulerian counterpart of the filter is the same filter without
advection (offline, `advection = nothing`), which uses the same ``M`` at fixed points; the Eulerian filter in
post-processing ([`compute_Eulerian_filter!`](@ref)) does not support a frequency scale.

## Mean velocities

The mean velocities are computed from the maps, and are the velocity of the mean position,
``\mathrm{d}\vb*{\varphi}^*/\mathrm{d}t``. The component-wise mean of the velocity, ``\vb*{u}^*``, is the same only
when the filter is time-invariant along the trajectory, i.e. its weights depend only on ``t - s``: then filtering
commutes with differentiation in time. That's the case for a uniform filter, or for a frequency scale that is
constant along the trajectory within the filter's time window. Otherwise, differentiating ``\vb*{\varphi}^*`` brings
out ``M`` at the current time, while ``\vb*{u}^*`` has ``M`` at each earlier time inside the average:

```math
\frac{\mathrm{d}\vb*{\varphi}^*}{\mathrm{d}t} = M(\vb*{\varphi}(t), t)\left(\frac{\vb*{u}}{M}\right)^{\!*},
\qquad
\frac{\mathrm{d}\vb*{\varphi}^*}{\mathrm{d}t} - \vb*{u}^* = \int_{-\infty}^{\infty} G\big(\tau(t)-\tau(s)\big)
\big[M(\vb*{\varphi}(t), t) - M(\vb*{\varphi}(s), s)\big]\,\vb*{u}(s)\,\mathrm{d}s ,
```

writing ``\vb*{\varphi}(s)`` and ``\vb*{u}(s)`` for the position and velocity of the particle at time ``s``. The
velocity of the mean position is consistent with the maps and with the regridding to the mean position. To get the
filtered velocity ``\vb*{u}^*`` instead, filter the velocity like a tracer: online, add it to `var_names_to_filter`
under another name (e.g. an auxiliary field kept equal to it); offline, save it under another name in the original
simulation.

## Online filtering

The online filter uses the forward part of these equations:
the filter clock advances at ``\mathrm{d}\tau/\mathrm{d}t=M``, and both the
tracer input and filter decay are multiplied by ``M``.

# Choosing online filters

## General form

The offline filter uses weight functions of the form
```math
G(t) = 
\begin{cases}
&\sum_{n=1}^{N/2} e^{-c_n t}(a_n\cos{d_n t} + b_n\sin{d_n t})\,, \hspace{1cm} &t > 0 \,, \\
& 0 &t \leq 0 \,,
\end{cases}
```
where ``a_n``, ``b_n``, ``c_n``, and ``d_n`` are real scalars, ``c_n > 0``, and ``N`` should be even. ``N`` is the number of exponentials that are summed to form the weight function, and should be even as the exponentials come in complex conjugate pairs to keep calculations real. These coefficients can be provided to [`OnlineFilterConfig`](@ref "OnlineFilterConfig") inside the `NamedTuple` `filter_params`.

```julia
filter_params = (a1 = 1, b1 = 1, c1 = 1, d1 = 1)
```
For the weight function to be normalised (so that the mean of a constant is the constant itself), these coefficients must be chosen such that

```math
\sum_{n=1}^{N/2} \frac{a_nc_n + b_n d_n}{c_n^2 +d_n^2} = 1\,.
```

Un-normalised filters can be used, but `map_to_mean` will be set to false as the maps no longer make sense. 

For ``N/2`` sets of coefficients, the weight function is composed of ``N`` exponentials, and ``N`` filtered tracers are needed to find the Lagrangian mean of each tracer. The number of equations that the filtering simulation solves is therefore linear in ``N``, so beware making ``N`` too large. 

## Exponential
For the special case of one (real) exponential, ``N`` can be set to 1 (this is the only exception to``N`` being even). The parameters ``a_1`` and ``c_1`` can then be provided:

```julia
filter_params = (a1 = 1, c1 = 1)
```
giving 
```math
G(t) = a_1 e^{-c_1 t}\,.
```

## Butterworth

Instead of providing the individual parameters in `filter_params`, the user can provide ``N`` (the filter order, which should be even or 1) and `freq_c` (the cut-off frequency) to use a filter with frequency response 

```math
\begin{equation}
    |\hat{G}(\omega)| = \frac{1}{\sqrt{1 + \left(\omega/\omega_c\right)^{2N}}}\,.
\end{equation} 
```

This is a Butterworth order-``N`` filter.


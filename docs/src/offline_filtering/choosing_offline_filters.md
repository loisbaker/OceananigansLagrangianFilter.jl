# Choosing offline filters

## General form

The offline filter uses weight functions of the form
```math
G(t) = \sum_{n=1}^{N/2} e^{-c_n|t|}(a_n\cos{d_n|t|} + b_n\sin{d_n|t|})\,,
```
where ``a_n``, ``b_n``, ``c_n``, and ``d_n`` are real scalars, ``c_n > 0``, and ``N`` should be even. ``N`` is the number of exponentials that are summed to form the weight function, and should be even as the exponentials come in complex conjugate pairs to keep calculations real. These coefficients can be provided to [`OfflineFilterConfig`](@ref "OfflineFilterConfig") inside the `NamedTuple` `filter_params`.

```julia
filter_params = (a1 = 0.5, b1 = 0.5, c1 = 1, d1 = 1)

```
For the weight function to be normalised (so that the mean of a constant is the constant itself), these coefficients must be chosen such that

```math
\sum_{n=1}^{N/2} \frac{a_nc_n + b_n d_n}{c_n^2 +d_n^2} = \frac{1}{2}\,.
```

Un-normalised filters can be used, but `map_to_mean` will be set to false as the maps no longer make sense. 

For ``N/2`` sets of coefficients, the weight function is composed of ``N`` exponentials, and ``N`` filtered tracers are needed to find the Lagrangian mean of each tracer. The number of equations that the filtering simulation solves is therefore linear in ``N``, so beware making ``N`` too large. 

## Exponential
For the special case of one (real) exponential, ``N`` can be set to 1 (this is the only exception to``N`` being even). The parameters ``a_1`` and ``c_1`` can then be provided:

```julia
filter_params = (a1 = 0.5, c1 = 1)
```
giving 
```math
G(t) = a_1 e^{-c_1 |t|}\,.
```

## Butterworth (squared)

Instead of providing the individual parameters in `filter_params`, the user can provide ``N`` (the filter order, which should be even or 1) and ``freq_c`` (the cut-off frequency) to use a filter with frequency response 

```math
\begin{equation}
    \hat{G}(\omega) = \frac{1}{1 + \left(\omega/\omega_c\right)^{2N}}\,.
\end{equation} 
```

This frequency response is the squared amplitude of that of the Butterworth order-``N`` filter.

The filter coefficients are set as
```math
\begin{align}
    a_n &= \frac{\omega_c}{N}\sin{\frac{\pi}{2N}(2n-1)}\,, \\
    b_n &= \frac{\omega_c}{N}\cos{\frac{\pi}{2N}(2n-1)}\,, \\
    c_n &= \omega_c\sin{\frac{\pi}{2N}(2n-1)}\,, \\
    d_n &= \omega_c\cos{\frac{\pi}{2N}(2n-1)}\,.
\end{align}
```




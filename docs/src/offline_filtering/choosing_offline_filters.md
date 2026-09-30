# Choosing offline filters

## General form

The offline filter uses weight functions of the form
```math
G(t) = \sum_{n=1}^{N/2} e^{-c_n|t|}(a_n\cos{d_n t} + b_n\sigma(t)\sin{d_n|t|})\,,
```
where ``a_n``, ``b_n``, ``c_n``, and ``d_n`` are real scalars, ``c_n > 0``, and ``N`` should be even. The sine terms are even in time by default, with ``\sigma(t) = 1``, so that ``G`` is even; they can instead be odd in time, with ``\sigma(t) = \mathrm{sgn}(t)`` (see `sine_parity` below). ``N`` is the number of exponentials that are summed to form the weight function, and should be even as the exponentials come in complex conjugate pairs to keep calculations real. These coefficients can be provided to [`OfflineFilterConfig`](@ref) inside the `NamedTuple` `filter_params`.

```julia
filter_params = (a1 = 0.5, b1 = 0.5, c1 = 1, d1 = 1)

```
`filter_params` can also contain two options for how the filter is applied and output, given as `Symbol`s:
- `sine_parity`: `:even` (default), for even sine terms, or `:odd`, for odd sine terms ``b_n \sin(d_n t)``.
- `outputs`: `:combined` (default), to output each filtered quantity as the sum of its terms (e.g. `b_Lagrangian_filtered`), or `:separate`, to output each term separately, named after its filtered tracer with a `_scaled` suffix (e.g. `b_C1_scaled` and `b_S1_scaled`, the ``a_1`` and ``b_1`` terms of the filtered ``b``, summed over the forward and backward passes). Separate outputs can't be regridded to the mean position.

For the weight function to be normalised (so that the mean of a constant is the constant itself, i.e. the gain at zero frequency is 1), these coefficients must be chosen such that

```math
\sum_{n=1}^{N/2} \frac{a_nc_n + b_n d_n}{c_n^2 +d_n^2} = \frac{1}{2}
```

for even sine terms. Odd sine terms integrate to zero, so for odd sine terms the ``b_n d_n`` terms are left out.

Un-normalised filters can be used (for example the spectral filter below), but `regrid_to_mean` will be set to false as the maps are no longer displacements from the mean position. 

To check a filter, [`get_weight_function`](@ref) and [`get_frequency_response`](@ref) compute its weight function and frequency response, for the whole filter or a single term (e.g. `term = "C1"`).

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

## Spectral filter

A spectral filter, [`set_offline_spectrum_filter_params`](@ref), extracts the signal in a narrow band around each of a set of frequencies ``\omega_n``, using an exponential window that lasts a time of about ``1/\alpha``:

```julia
filter_params = set_offline_spectrum_filter_params(alpha = 1e-5, freqs = [1e-4, 2e-4], normalisation = :unit_gain)
```

The weight functions at each frequency are
```math
\begin{align}
    G_{Cn}(t) &= A_n e^{-\alpha|t|}\cos{\omega_n t}\,, \\
    G_{Sn}(t) &= A_n e^{-\alpha|t|}\sin{\omega_n t}\,,
\end{align}
```
so that ``a_n = b_n = A_n``, ``c_n = \alpha`` and ``d_n = \omega_n``, with odd sine terms (`sine_parity = :odd`). Each term is output separately (`outputs = :separate`), e.g. `b_C1_scaled` and `b_S1_scaled` for a filtered variable `b` at ``\omega_1``.

There are two normalisations:
- `:spectral` (default): ``A_n = \sqrt{\alpha}``, so that the window has unit energy. Then the sum of the squares of the two outputs at ``\omega_n`` (e.g. `b_C1_scaled^2 + b_S1_scaled^2`) estimates the (two-sided) power spectral density of the signal at ``\omega_n``.
- `:unit_gain`: ``A_n = \alpha(\alpha^2 + 4\omega_n^2)/(2\alpha^2 + 4\omega_n^2)``, so that the cosine term (e.g. `b_C1_scaled`) passes a signal at exactly ``\omega_n`` through unchanged. Its frequency response is
```math
\hat{G}_{Cn}(\omega) = A_n\left(\frac{\alpha}{\alpha^2 + (\omega - \omega_n)^2} + \frac{\alpha}{\alpha^2 + (\omega + \omega_n)^2}\right)\,,
```
which is real (the weight function is even), and ``A_n`` is chosen so that ``\hat{G}_{Cn}(\omega_n) = 1``.

The Eulerian filter for comparison (`compute_Eulerian_filter = true`) also outputs the terms separately, e.g. `b_Eulerian_filtered_C1_scaled`.

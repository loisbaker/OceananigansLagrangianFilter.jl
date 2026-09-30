# Choosing online filters

## General form

The online filter uses weight functions of the form
```math
G(t) = 
\begin{cases}
&\sum_{n=1}^{N/2} e^{-c_n t}(a_n\cos{d_n t} + b_n\sin{d_n t})\,, \hspace{1cm} &t > 0 \,, \\
& 0 &t \leq 0 \,,
\end{cases}
```
where ``a_n``, ``b_n``, ``c_n``, and ``d_n`` are real scalars, ``c_n > 0``, and ``N`` should be even. ``N`` is the number of exponentials that are summed to form the weight function, and should be even as the exponentials come in complex conjugate pairs to keep calculations real. These coefficients can be provided to [`OnlineFilterConfig`](@ref) inside the `NamedTuple` `filter_params`.

```julia
filter_params = (a1 = 1, b1 = 1, c1 = 1, d1 = 1)
```
For the weight function to be normalised (so that the mean of a constant is the constant itself), these coefficients must be chosen such that

```math
\sum_{n=1}^{N/2} \frac{a_nc_n + b_n d_n}{c_n^2 +d_n^2} = 1\,.
```

Un-normalised filters can be used (for example the spectral filter below), but `regrid_to_mean` will be set to false as the maps are no longer displacements from the mean position. 

`filter_params` can also contain the option `outputs`: `:combined` (default), to output each filtered quantity as the sum of its terms (e.g. `b_Lagrangian_filtered`), or `:separate`, to output each term separately, named after its filtered tracer with a `_scaled` suffix (e.g. `b_C1_scaled` and `b_S1_scaled`, the ``a_1`` and ``b_1`` terms of the filtered ``b``). Separate outputs can't be regridded to the mean position. (The option `sine_parity` of the offline filter has no effect here, since the online weight function is only used for ``t > 0``.)

To check a filter, [`get_weight_function`](@ref) and [`get_frequency_response`](@ref) (with `direction = :forward`) compute its weight function and frequency response.

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

## Spectral filter

A spectral filter, [`set_online_spectrum_filter_params`](@ref), is the causal version of the offline spectral filter (see [Choosing offline filters](@ref)). It extracts the signal in a narrow band around each of a set of frequencies ``\omega_n``, using an exponential window over the past that lasts a time of about ``1/\alpha``:

```julia
filter_params = set_online_spectrum_filter_params(alpha = 1e-5, freqs = [1e-4, 2e-4], normalisation = :unit_gain)
```

The weight functions at each frequency are, for ``t > 0``,
```math
\begin{align}
    G_{Cn}(t) &= A_n e^{-\alpha t}\cos{\omega_n t}\,, \\
    G_{Sn}(t) &= A_n e^{-\alpha t}\sin{\omega_n t}\,,
\end{align}
```
so that ``a_n = b_n = A_n``, ``c_n = \alpha`` and ``d_n = \omega_n``. Each term is output separately (`outputs = :separate`), e.g. `b_C1_scaled` and `b_S1_scaled` for a filtered variable `b` at ``\omega_1``. Because the window only uses the past, the outputs respond to changes in the signal with a delay of about ``1/\alpha``.

There are two normalisations:
- `:spectral` (default): ``A_n = \sqrt{2\alpha}``, so that the (one-sided) window has unit energy. Then, as for the offline filter, the sum of the squares of the two outputs at ``\omega_n`` (e.g. `b_C1_scaled^2 + b_S1_scaled^2`) estimates the (two-sided) power spectral density of the signal at ``\omega_n``.
- `:unit_gain`: ``A_n = \alpha\sqrt{\alpha^2 + 4\omega_n^2}/\sqrt{\alpha^2 + \omega_n^2}``, so that the cosine term (e.g. `b_C1_scaled`) passes a signal at exactly ``\omega_n`` with its amplitude preserved. Its frequency response is
```math
\hat{G}_{Cn}(\omega) = \frac{A_n}{2}\left(\frac{1}{\alpha + i(\omega - \omega_n)} + \frac{1}{\alpha + i(\omega + \omega_n)}\right)\,,
```
and ``A_n`` is chosen so that ``|\hat{G}_{Cn}(\omega_n)| = 1``. Unlike the offline filter, ``\hat{G}_{Cn}(\omega_n)`` is complex, since the weight function is one-sided, so the phase of the signal is shifted slightly, by about ``\alpha/(2\omega_n)`` radians.

The Eulerian filter for comparison ([`compute_Eulerian_filter!`](@ref)) also outputs the terms separately, e.g. `b_Eulerian_filtered_C1_scaled`.

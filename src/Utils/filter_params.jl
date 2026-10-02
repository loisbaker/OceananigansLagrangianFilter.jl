# Filter parameters: the constructors for the filter coefficients, helpers for the terms of a filter, and its weight
# function and frequency response

"""
    set_offline_BW2_filter_params(; N::Int=1, freq_c::Real=1)

Calculates the coefficients for a filter that has a frequency response given by a
Butterworth filter with order `N` and cutoff frequency `freq_c`, squared. 

Uses `N` exponentials and `N/2` sets of coefficients (a,b,c,d). N should therefore be even,
since exponentials come in pairs to ensure a real-valued filter.

However, the special case N=1 is allowed, which gives a single (real) exponential filter.

Frequency response: `Ghat(omega) = 1 / (1 + (omega / freq_c)^(2*N))`
Real filter shape: `G(t) = sum_{i=1}^{N/2} exp(-c_i*abs(t))*(a_i*cos(d_i * abs(t)) + b_i*sin(d_i * abs(t)))`

This function supports two types of filters:
* A **single exponential filter** when `N=1`. This is a special case that
  generates two coefficients instead of 4. The unidirectional filter is a single exponential,
  and `N_coeffs = 0.5`. Only `a1` and `c1` are returned.
* A **Butterworth squared filter** for `N>1`. This generates `N/2` sets of
  coefficients (a,b,c,d), representing a filter of order `N`. The
  coefficients are computed based on the filter's order and cutoff frequency.

Arguments
=========
- `N`: The order parameter for the filter. `N=1` for a single exponential.
  For `N>1`, the filter's order is `N`. Must be a non-negative even integer.
- `freq_c`: The cutoff frequency of the filter. Must be a real number.

Returns
=======
- A `NamedTuple` containing the filter coefficients, `N_coeffs` (the number
  of coefficient pairs), and the output options `outputs = :combined` and
  `sine_parity = :even` (the defaults, see `filter_output_options`).
"""
function set_offline_BW2_filter_params(;N::Int=1,freq_c::Real=1) 
    if N == 1
        # Special case single exponential
        N_coeffs = 0.5
    elseif N/2 != floor(N/2) || N < 0 
        error("N must be a non-negative even integer, or 1 for a single exponential.")
    else
        N_coeffs = Int(N/2)
    end

    freq_c = abs(freq_c) # Ensure freq_c is positive
    if N_coeffs == 0.5 # special case N=1, single exponential only has a cosine component
        a1 = freq_c/2 
        c1 = freq_c
        filter_params = (; a1 = a1, c1 = c1, N_coeffs = N_coeffs, outputs = :combined, sine_parity = :even)
        return filter_params

    else
        filter_params = NamedTuple()
        for i in 1:N_coeffs
            
            a = (freq_c/N)*sin(pi/2/N*(2*i-1))
            b = (freq_c/N)*cos(pi/2/N*(2*i-1))
            c = freq_c*sin(pi/2/N*(2*i-1))
            d = freq_c*cos(pi/2/N*(2*i-1))

            temp_params = NamedTuple{(Symbol("a$i"), Symbol("b$i"),Symbol("c$i"),Symbol("d$i"))}([a,b,c,d])
            filter_params = merge(filter_params,temp_params)
        end

        return merge(filter_params, (; N_coeffs = N_coeffs, outputs = :combined, sine_parity = :even))
    end
end

"""
    set_online_BW_filter_params(; N::Int=1, freq_c::Real=1)

Calculates the coefficients for a filter that has a frequency response given by a
Butterworth filter with order `N` and cutoff frequency `freq_c`. Note that the frequency
response is not squared, like in the offline forward-backward filter, and the frequency response
is not real-valued, implying a nonlinear phase shift. 

Uses `N` exponentials and `N/2` sets of coefficients (a,b,c,d). N should therefore be even,
since exponentials come in pairs to ensure a real-valued filter.

However, the special case N=1 is allowed, which gives a single (real) exponential filter.

Frequency response: `abs(Ghat(omega)) = 1 / sqrt(1 + (omega / freq_c)^(2*N))`
Real filter shape: `G(t) = sum_{i=1}^{N/2} exp(-c_i*t)* (a_i*cos(d_i * t) + b_i*sin(d_i * t))` 
for t>=0, and 0 for t<0.

This function supports two types of filters:
* A **single exponential filter** when `N=1`. This is a special case that
  generates two coefficients instead of 4. The unidirectional filter is a single exponential,
  and `N_coeffs = 0.5`. Only `a1` and `c1` are returned.
* A **Butterworth filter** for `N>1`. This generates `N/2` sets of
  coefficients (a,b,c,d), representing a filter of order `N`. The
  coefficients are computed based on the filter's order and cutoff frequency.

Arguments
=========
- `N`: The order parameter for the filter. `N=1` for a single exponential.
  For `N>1`, the filter's order is `N`. Must be a non-negative even integer.
- `freq_c`: The cutoff frequency of the filter. Must be a real number.

Returns
=======
- A `NamedTuple` containing the filter coefficients, `N_coeffs` (the number
  of coefficient pairs), and the output option `outputs = :combined` (the default,
  see `filter_output_options`).
"""
function set_online_BW_filter_params(;N::Int=1,freq_c::Real=1) 
    if N == 1
        # Special case single exponential
        N_coeffs = 0.5
    elseif N/2 != floor(N/2) || N < 0 
        error("N must be a non-negative even integer, or 1 for a single exponential.")
    else
        N_coeffs = Int(N/2)
    end

    freq_c = abs(freq_c) # Ensure freq_c is positive
    if N_coeffs == 0.5 # special case N=1, single exponential only has a cosine component
        a1 = freq_c
        c1 = freq_c
        filter_params = (; a1 = a1, c1 = c1, N_coeffs = N_coeffs, outputs = :combined)
        return filter_params

    else
        filter_params = NamedTuple()
        for i in 1:N_coeffs
            
            c = freq_c*sin(pi/2/N*(2*i-1))
            d = -freq_c*cos(pi/2/N*(2*i-1))
            ri = exp(pi*im*(2*i-1)/ (2*N)) # ith of the 2N roots of -1
            Ai = 1
            for k in 1:N
                if k != i
                    rk = exp(pi*im*(2*k-1)/ (2*N))
                    Ai *= 1/(ri - rk)
                end
            end
            b = 2*freq_c*real(Ai*exp(-im*pi*N/2))
            a = -2*freq_c*imag(Ai*exp(-im*pi*N/2))

            temp_params = NamedTuple{(Symbol("a$i"), Symbol("b$i"),Symbol("c$i"),Symbol("d$i"))}([a,b,c,d])
            filter_params = merge(filter_params,temp_params)
        end

        return merge(filter_params, (; N_coeffs = N_coeffs, outputs = :combined))
    end
end

"""
    set_offline_spectrum_filter_params(; alpha::Real, freqs::AbstractVector, normalisation::Symbol = :spectral)

Coefficients for an exponentially-windowed spectral filter, which extracts the signal in a narrow band
around each frequency in `freqs`. The window is `w(t) = exp(-alpha*|t|)`, giving the in-phase (`C`) and
quadrature (`S`) weight functions at each frequency `omega_n`:

    C(t) = A_n exp(-alpha*|t|) cos(omega_n*t)
    S(t) = A_n exp(-alpha*|t|) sin(omega_n*t)

so that `a_n = b_n = A_n`, `c_n = alpha` and `d_n = omega_n`. The `S` weight functions are odd in `t`
(`sine_parity = :odd`), and each term is output separately (`outputs = :separate`).

`alpha` sets the frequency resolution: the window lasts about `1/alpha`, so frequencies closer than
about `alpha` are not resolved.

Two normalisations are available:
- `:spectral` (default): `A_n = sqrt(alpha)`, so that the window has unit energy, `∫ (A_n w)^2 dt = 1`.
  Then `C^2 + S^2` estimates the (two-sided) power spectral density of the signal at `omega_n`,
  normalised so that the variance is `∫ S(omega) domega / 2π`.
- `:unit_gain`: `A_n = alpha (alpha^2 + 4omega_n^2) / (2alpha^2 + 4omega_n^2)`, so that `C` passes a signal at
  exactly `omega_n` through unchanged: `C` is the band-passed signal, and `S` is its quadrature (a quarter
  period later). Then `sqrt(C^2 + S^2)` is its envelope and `atan(S, C)` its instantaneous phase (`omega_n*t`
  plus the phase of the signal). `S` has gain `4omega_n^2 / (2alpha^2 + 4omega_n^2)`, so the envelope and phase
  are accurate to about `alpha^2 / (2omega_n^2)`. At `omega_n = 0` this reduces to `alpha/2`, the single
  exponential of `set_offline_BW2_filter_params(N = 1)`.

The `S` terms can be left out of the outputs by setting their coefficients `b_n` to zero. See
[`set_online_spectrum_filter_params`](@ref) for the online filter.

Arguments
=========
- `alpha`: Decay rate of the exponential window. Must be positive.
- `freqs`: Frequencies (radians per unit time) at which to extract the signal.
- `normalisation`: `:spectral` (default) or `:unit_gain`.

Returns
=======
- A `NamedTuple` of coefficients, `N_coeffs`, `outputs = :separate` and `sine_parity = :odd`.
"""
function set_offline_spectrum_filter_params(; alpha::Real, freqs::AbstractVector, normalisation::Symbol = :spectral)
    amplitude(omega) = normalisation === :spectral ? sqrt(alpha) :
                       alpha * (alpha^2 + 4omega^2) / (2alpha^2 + 4omega^2) # Unit gain at omega
    filter_params = spectrum_filter_coefficients(amplitude; alpha, freqs, normalisation)
    return merge(filter_params, (; sine_parity = :odd))
end

"""
    set_online_spectrum_filter_params(; alpha::Real, freqs::AbstractVector, normalisation::Symbol = :spectral)

Coefficients for an exponentially-windowed spectral filter for the online filter: the causal version of
[`set_offline_spectrum_filter_params`](@ref), which only uses the past. The window is `w(t) = exp(-alpha*t)`
for `t > 0` (the time before the present), giving the in-phase (`C`) and quadrature (`S`) weight functions
at each frequency `omega_n`:

    C(t) = A_n exp(-alpha*t) cos(omega_n*t)
    S(t) = A_n exp(-alpha*t) sin(omega_n*t)

so that `a_n = b_n = A_n`, `c_n = alpha` and `d_n = omega_n`, and each term is output separately
(`outputs = :separate`).

`alpha` sets the frequency resolution: the window lasts about `1/alpha`, so frequencies closer than
about `alpha` are not resolved. Because the window only uses the past, changes in the amplitude of the
signal show up in the envelope `sqrt(C^2 + S^2)` after a delay of about `1/alpha` (so `t - 1/alpha` is the time
for the envelope or power, rather than the time computed by [`compute_time_shift!`](@ref)).

Two normalisations are available:
- `:spectral` (default): `A_n = sqrt(2alpha)`, so that the window has unit energy, `∫ (A_n w)^2 dt = 1`.
  Then, as for the offline filter, `C^2 + S^2` estimates the (two-sided) power spectral density of the
  signal at `omega_n`, normalised so that the variance is `∫ S(omega) domega / 2π`.
- `:unit_gain`: `A_n = alpha sqrt(alpha^2 + 4omega_n^2) / sqrt(alpha^2 + omega_n^2)`, so that `C` passes a
  signal at exactly `omega_n` with its amplitude preserved, and `S` gives its quadrature with gain
  `omega_n / sqrt(alpha^2 + omega_n^2)`. Then `sqrt(C^2 + S^2)` is its envelope and `atan(S, C)` its
  instantaneous phase (`omega_n*t` plus the phase of the signal). Unlike the offline filter, the one-sided
  window shifts the phases slightly: `C` lags the signal and `S` leads its quadrature, each by about
  `alpha / (2omega_n)` radians, so the envelope and phase have a ripple of relative size about
  `alpha / (2omega_n)` at twice the frequency (for a ripple below 1%, use `alpha ≲ 0.02 omega_n`). At
  `omega_n = 0` this reduces to `alpha`, the single exponential of `set_online_BW_filter_params(N = 1)`.

The `S` terms can be left out of the outputs by setting their coefficients `b_n` to zero.

Arguments
=========
- `alpha`: Decay rate of the exponential window. Must be positive.
- `freqs`: Frequencies (radians per unit time) at which to extract the signal.
- `normalisation`: `:spectral` (default) or `:unit_gain`.

Returns
=======
- A `NamedTuple` of coefficients, `N_coeffs` and `outputs = :separate`.
"""
function set_online_spectrum_filter_params(; alpha::Real, freqs::AbstractVector, normalisation::Symbol = :spectral)
    amplitude(omega) = normalisation === :spectral ? sqrt(2alpha) :
                       alpha * sqrt(alpha^2 + 4omega^2) / sqrt(alpha^2 + omega^2) # Unit gain at omega
    return spectrum_filter_coefficients(amplitude; alpha, freqs, normalisation)
end

# Coefficients of a spectral filter with an exponential window decaying at rate alpha, with a_n = b_n = amplitude(omega_n),
# c_n = alpha and d_n = omega_n at each frequency omega_n in freqs, whose terms are output separately
function spectrum_filter_coefficients(amplitude::Function; alpha::Real, freqs::AbstractVector, normalisation::Symbol)
    alpha > 0 || error("alpha must be positive.")
    length(freqs) > 0 || error("freqs must contain at least one frequency.")
    normalisation in (:spectral, :unit_gain) || error("normalisation must be :spectral or :unit_gain, got $(repr(normalisation))")

    filter_params = NamedTuple()
    for (n, omega) in enumerate(freqs)
        A = amplitude(omega)
        coefficients = NamedTuple{(Symbol("a$n"), Symbol("b$n"), Symbol("c$n"), Symbol("d$n"))}((A, A, alpha, omega))
        filter_params = merge(filter_params, coefficients)
    end
    return merge(filter_params, (; N_coeffs = length(freqs), outputs = :separate))
end

# The number of filter coefficients (a1, b1, c1, d1, a2, ...) in filter_params, not counting its other fields
number_of_coefficient_fields(filter_params::NamedTuple) = count(k -> occursin(r"^[abcd][0-9]+$", String(k)), keys(filter_params))

"""
    filter_output_options(filter_params::NamedTuple)

Return the output options `(outputs, sine_parity)` of a filter, with their defaults if not given:

- `outputs`: `:combined` (default), to output each filtered quantity as a single field, or `:separate`, to
  output each (cosine or sine) term of the filter separately, named after its tracer with a `_scaled` suffix.
- `sine_parity`: `:even` (default), if the sine terms of the weight function are `sin(d|t|)`, or `:odd`, if they
  are `sin(dt)`. Odd sine terms change sign in the backward pass of the offline filter.
"""
function filter_output_options(filter_params::NamedTuple)
    outputs     = get(filter_params, :outputs, :combined)
    sine_parity = get(filter_params, :sine_parity, :even)
    outputs in (:combined, :separate) || error("filter_params.outputs must be :combined or :separate, got $(repr(outputs))")
    sine_parity in (:even, :odd)      || error("filter_params.sine_parity must be :even or :odd, got $(repr(sine_parity))")
    return outputs, sine_parity
end

# The terms (i, has_cosine, has_sine) of a filter that are output separately: terms with a zero coefficient are skipped
filter_terms(filter_params) = filter_params.N_coeffs == 0.5 ? [(1, filter_params.a1 != 0, false)] :
    [(i, getproperty(filter_params, Symbol("a$i")) != 0, getproperty(filter_params, Symbol("b$i")) != 0) for i in 1:filter_params.N_coeffs]

# The coefficients of each term of a filter that is output separately, as a single-term filter, with the suffix of
# its output names (see filtered_output_names)
function separate_term_params(filter_params)
    _, sine_parity = filter_output_options(filter_params)
    terms = Tuple{String, NamedTuple}[]
    for (i, has_cosine, has_sine) in filter_terms(filter_params)
        a, c = getproperty(filter_params, Symbol("a$i")), getproperty(filter_params, Symbol("c$i"))
        if filter_params.N_coeffs == 0.5
            push!(terms, ("_C1_scaled", (; a1 = a, c1 = c, N_coeffs = 0.5)))
        else
            b, d = getproperty(filter_params, Symbol("b$i")), getproperty(filter_params, Symbol("d$i"))
            has_cosine && push!(terms, ("_C$(i)_scaled", (; a1 = a, b1 = zero(b), c1 = c, d1 = d, N_coeffs = 1, sine_parity)))
            has_sine   && push!(terms, ("_S$(i)_scaled", (; a1 = zero(a), b1 = b, c1 = c, d1 = d, N_coeffs = 1, sine_parity)))
        end
    end
    return terms
end

# The coefficients of one term of a filter, labelled as in the names of its separate outputs (e.g. term = "C1" for
# b_C1_scaled, or "S2"), or of the whole filter if term = nothing
function select_term(filter_params::NamedTuple, term)
    isnothing(term) && return filter_params
    terms = Dict(suffix => params for (suffix, params) in separate_term_params(filter_params))
    key = "_" * String(term) * "_scaled"
    haskey(terms, key) || error("term must be one of $(join(sort([k[2:end-7] for k in keys(terms)]), ", ")) (or nothing, " *
                                "for the whole filter), got $(repr(term))")
    return terms[key]
end

"""
    get_weight_function(; t::AbstractArray, tref::Real, filter_params::NamedTuple, direction::Symbol = :both, term = nothing)

Computes the weight function of the filter, which determines how much each point in the
timeseries `t` contributes to the filtered value at a reference time `tref`. For the offline
filter this is

``G(t) = \\sum_i e^{-c_i |\\tau|} \\left[a_i \\cos(d_i \\tau) + b_i \\sin(d_i |\\tau|)\\right], \\quad \\tau = t_{ref} - t,``

or with sine terms ``b_i \\sin(d_i \\tau)`` if `filter_params.sine_parity = :odd`.

Keyword arguments
=================
- `t`: A collection of time points in the timeseries.
- `tref`: The reference time at which the filter is being evaluated.
- `filter_params`: A `NamedTuple` containing the coefficients (`a`, `b`, `c`,
  `d`) and the number of coefficient pairs (`N_coeffs`).
- `direction`: `:both` (default) for the offline filter, `:forward` to keep only times before `tref`
  (as for the online filter, or the forward pass of the offline filter), or `:backward` to keep only
  times after `tref` (the backward pass of the offline filter).
- `term`: `nothing` (default) for the whole filter, or a single term, labelled as in the names of the separate
  outputs (see [`create_output_fields`](@ref)), e.g. `"C1"` (for `b_C1_scaled`) or `"S2"`.

Returns
=======
- A vector of weights `G`, with the same dimensions as `t`, representing the
  value of the weight function at each time point relative to `tref`.
"""
function get_weight_function(; t::AbstractArray, tref::Real, filter_params::NamedTuple, direction::Symbol = :both, term = nothing)
    filter_params = select_term(filter_params, term)
    direction in (:both, :forward, :backward) || error("direction must be :both, :forward or :backward, got $(repr(direction))")
    _, sine_parity = filter_output_options(filter_params)
    τ = tref .- t # Time lag, positive for times before tref

    G = 0*t
    N_coeffs = filter_params.N_coeffs
    if N_coeffs == 0.5
        a1 = filter_params.a1
        c1 = filter_params.c1
        G .= a1.*exp.(-c1.*abs.(τ))
    else
        for i in 1:N_coeffs

            a = getproperty(filter_params, Symbol("a$i"))
            b = getproperty(filter_params, Symbol("b$i"))
            c = getproperty(filter_params, Symbol("c$i"))
            d = getproperty(filter_params, Symbol("d$i"))

            sine = sine_parity === :odd ? sin.(d.*τ) : sin.(d.*abs.(τ))
            G += (a.*cos.(d.*abs.(τ)) .+ b.*sine).*exp.(-c.*abs.(τ))
        end
    end
    if direction === :forward
        G[t .> tref] .= 0
    elseif direction === :backward
        G[t .< tref] .= 0
    end
    return G
end

"""
    get_frequency_response(; freq::AbstractArray, filter_params::NamedTuple, direction::Symbol = :both, term = nothing)

Calculates the frequency response ``\\hat{G}(\\omega)`` of the filter, i.e. the Fourier transform of the
weight function (see [`get_weight_function`](@ref)), so that filtering ``e^{i\\omega t}`` gives
``\\hat{G}(\\omega) e^{i\\omega t}``.

Keyword arguments
=================
- `freq`: A vector of frequencies ``\\omega`` (in radians per unit time).
- `filter_params`: A `NamedTuple` containing the filter coefficients (`a`, `b`,
  `c`, `d`) and the number of coefficient pairs (`N_coeffs`).
- `direction`: `:both` (default) for the offline filter, `:forward` for the online filter (or the
  forward pass of the offline filter), or `:backward` for the backward pass of the offline filter.
- `term`: `nothing` (default) for the whole filter, or a single term, labelled as in the names of the separate
  outputs (see [`create_output_fields`](@ref)), e.g. `"C1"` (for `b_C1_scaled`) or `"S2"`.

Returns
=======
- A vector `Ghat` of the frequency response at each frequency in `freq`. This is real for the offline
  filter (`direction = :both`) if the weight function is even in time (`sine_parity = :even`, the
  default), and complex otherwise.
"""
function get_frequency_response(; freq::AbstractArray, filter_params::NamedTuple, direction::Symbol = :both, term = nothing)
    filter_params = select_term(filter_params, term)
    direction in (:both, :forward, :backward) || error("direction must be :both, :forward or :backward, got $(repr(direction))")
    _, sine_parity = filter_output_options(filter_params)
    σ = sine_parity === :odd ? -1 : 1 # Sign of the sine terms for times after tref

    # The weight function for times after tref is that for times before tref reversed in time, with the sine terms multiplied by σ
    Ghat = zeros(Complex{float(eltype(freq))}, size(freq))
    direction in (:both, :forward)  && (Ghat .+= one_sided_frequency_response(freq, filter_params, 1))
    direction in (:both, :backward) && (Ghat .+= one_sided_frequency_response(-freq, filter_params, σ))

    # The response of the whole weight function is real if it is even in time
    return (direction === :both && sine_parity === :even) ? real(Ghat) : Ghat
end

# The frequency response of the weight function for times before tref, with the sine terms multiplied by σ
function one_sided_frequency_response(freq, filter_params, σ)
    if filter_params.N_coeffs == 0.5
        return filter_params.a1 ./ (filter_params.c1 .+ im .* freq)
    end
    Ghat = zeros(Complex{float(eltype(freq))}, size(freq))
    for i in 1:filter_params.N_coeffs
        a = getproperty(filter_params, Symbol("a$i"))
        b = σ * getproperty(filter_params, Symbol("b$i"))
        c = getproperty(filter_params, Symbol("c$i"))
        d = getproperty(filter_params, Symbol("d$i"))
        Ghat .+= 0.5 .* ((a + im*b) ./ (c .+ im .* (d .+ freq)) .+ (a - im*b) ./ (c .+ im .* (freq .- d)))
    end
    return Ghat
end

# The gain of the filter at zero frequency (the integral of its weight function), which is 1 for a normalised filter
zero_frequency_gain(filter_params::NamedTuple, direction::Symbol) =
    real(get_frequency_response(freq = [0.0]; filter_params, direction)[1])

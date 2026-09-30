@testset "OfflineFilterConfig boundary relaxation validation" begin

    good_mask(x, z, p) = 1.0 # one argument per non-Flat dimension, plus mask_params

    # boundary_relaxation = true requires relax_timescale to be set
    @test_throws ErrorException OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                    var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                    N = 1, freq_c = 1e-4,
                                    boundary_relaxation = true,
                                    mask_func = good_mask, mask_params = (; a = 1))

    # boundary_relaxation = true requires a mask_func to be set
    @test_throws ErrorException OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                    var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                    N = 1, freq_c = 1e-4,
                                    boundary_relaxation = true,
                                    relax_timescale = 1hour, mask_params = (; a = 1))

    # mask_func must take one argument per non-Flat grid dimension (plus mask_params)
    bad_mask(x, p) = 1.0 # missing the z argument
    @test_throws ErrorException OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                    var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                    N = 1, freq_c = 1e-4,
                                    boundary_relaxation = true, relax_timescale = 1hour,
                                    mask_func = bad_mask, mask_params = (; a = 1))

    # mask_params = nothing is only a warning, not an error, since mask_func might not need parameters
    config = @test_logs (:warn,) match_mode=:any OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                    var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                    N = 1, freq_c = 1e-4,
                                    boundary_relaxation = true, relax_timescale = 1hour,
                                    mask_func = good_mask, mask_params = nothing)
    @test config.boundary_relaxation
    @test config.relax_timescale == 1hour

    # A fully-specified relaxation configuration constructs without error or warning
    config2 = OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                    var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                    N = 1, freq_c = 1e-4,
                                    boundary_relaxation = true, relax_timescale = 1hour,
                                    mask_func = good_mask, mask_params = (; a = 1))
    @test config2.mask_func === good_mask
    @test config2.mask_params == (; a = 1)
end

@testset "OnlineFilterConfig boundary relaxation validation" begin

    grid = RectilinearGrid(CPU(), size = (4, 4), x = (-1, 1), z = (-1, 0),
                            topology = (Periodic, Flat, Bounded))

    good_mask(x, z, p) = 1.0

    @test_throws ErrorException OnlineFilterConfig(grid = grid,
                                    var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                    N = 1, freq_c = 1e-4,
                                    boundary_relaxation = true,
                                    mask_func = good_mask, mask_params = (; a = 1))

    @test_throws ErrorException OnlineFilterConfig(grid = grid,
                                    var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                    N = 1, freq_c = 1e-4,
                                    boundary_relaxation = true,
                                    relax_timescale = 1hour, mask_params = (; a = 1))

    bad_mask(x, p) = 1.0
    @test_throws ErrorException OnlineFilterConfig(grid = grid,
                                    var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                    N = 1, freq_c = 1e-4,
                                    boundary_relaxation = true, relax_timescale = 1hour,
                                    mask_func = bad_mask, mask_params = (; a = 1))

    config = @test_logs (:warn,) match_mode=:any OnlineFilterConfig(grid = grid,
                                    var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                    N = 1, freq_c = 1e-4,
                                    boundary_relaxation = true, relax_timescale = 1hour,
                                    mask_func = good_mask, mask_params = nothing)
    @test config.boundary_relaxation
    @test config.relax_timescale == 1hour

    config2 = OnlineFilterConfig(grid = grid,
                                    var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                    N = 1, freq_c = 1e-4,
                                    boundary_relaxation = true, relax_timescale = 1hour,
                                    mask_func = good_mask, mask_params = (; a = 1))
    @test config2.mask_func === good_mask
    @test config2.mask_params == (; a = 1)
end

@testset "Map options (compute_maps, regrid_to_mean)" begin

    base = (; original_data_filename = "data/reference_sim.jld2", var_names_to_filter = ("b",),
              velocity_names = ("u", "w"), N = 1, freq_c = 1e-4)

    # By default, the maps are computed and the filtered fields are regridded to the mean position
    config = OfflineFilterConfig(; base...)
    @test config.compute_maps && config.regrid_to_mean && config.compute_mean_velocities

    # The maps can be output without regridding
    config = OfflineFilterConfig(; base..., compute_maps = true, regrid_to_mean = false)
    @test config.compute_maps && !config.regrid_to_mean

    # Regridding requires the maps, so they are turned on with a warning
    config = @test_logs (:warn, r"requires the maps") match_mode=:any OfflineFilterConfig(; base..., compute_maps = false)
    @test config.compute_maps && config.regrid_to_mean

    # For the Eulerian filter the maps and mean velocities can still be computed, but not regridded
    config = @test_logs (:warn, r"Eulerian") match_mode=:any OfflineFilterConfig(; base..., advection = nothing)
    @test config.compute_maps && !config.regrid_to_mean && config.compute_mean_velocities

    # Likewise for un-normalised filter coefficients
    config = @test_logs (:warn, r"not normalised") match_mode=:any OfflineFilterConfig(; base..., N = nothing, freq_c = nothing,
                                                                                         filter_params = (a1 = 1.0, c1 = 1.0))
    @test config.compute_maps && !config.regrid_to_mean

    # The removed option map_to_mean gives an error explaining what to use instead
    @test_throws r"compute_maps" OfflineFilterConfig(; base..., map_to_mean = true)
end

@testset "Filter output options (outputs, sine_parity)" begin

    base = (; original_data_filename = "data/reference_sim.jld2", var_names_to_filter = ("b",), velocity_names = ("u", "w"))
    butterworth = set_offline_BW2_filter_params(N = 2, freq_c = 1e-4)

    # The helpers set the output options explicitly (these are also the defaults)
    @test butterworth.outputs === :combined && butterworth.sine_parity === :even
    @test set_online_BW_filter_params(N = 2, freq_c = 1e-4).outputs === :combined

    # Only the allowed values can be used
    @test_throws r"outputs must be" OfflineFilterConfig(; base..., filter_params = merge(butterworth, (outputs = :both,)))
    @test_throws r"sine_parity must be" OfflineFilterConfig(; base..., filter_params = merge(butterworth, (sine_parity = :none,)))

    # A single exponential filter has no sine terms
    @test_throws r"no sine terms|has none" OfflineFilterConfig(; base..., filter_params = merge(set_offline_BW2_filter_params(N = 1, freq_c = 1e-4), (sine_parity = :odd,)))

    # N_coeffs can be inferred from the coefficients, even if filter_params also has the output options
    without_N_coeffs = Base.structdiff(butterworth, (; N_coeffs = nothing))
    @test OfflineFilterConfig(; base..., filter_params = without_N_coeffs).filter_params.N_coeffs == 1

    # Separate outputs can't be regridded
    separate = merge(butterworth, (outputs = :separate,))
    config = @test_logs (:warn, r"no combined fields") match_mode=:any OfflineFilterConfig(; base..., filter_params = separate)
    @test config.compute_maps && !config.regrid_to_mean

    # Odd sine terms integrate to zero, so they don't count towards the normalisation: this term has
    # gain 2ac/(c^2 + d^2) = 1 with odd sine terms, but 2(ac + bd)/(c^2 + d^2) = 2 with even ones
    single_term = (; a1 = 1.0, b1 = 1.0, c1 = 1.0, d1 = 1.0, N_coeffs = 1, outputs = :combined)
    @test_logs (:warn, r"gain at zero frequency is 2.0") match_mode=:any OfflineFilterConfig(; base...,
                    filter_params = merge(single_term, (; sine_parity = :even)))
    config = OfflineFilterConfig(; base..., filter_params = merge(single_term, (; sine_parity = :odd)))
    @test config.regrid_to_mean # not switched off, so it was treated as normalised

    # The online filter only uses the weight function for t > 0, so the sine parity has no effect
    grid = RectilinearGrid(size = (4, 4), x = (0, 1), z = (-1, 0), topology = (Periodic, Flat, Bounded))
    online_params = merge(set_online_BW_filter_params(N = 2, freq_c = 1e-4), (sine_parity = :odd,))
    @test_logs (:warn, r"no effect for the online filter") match_mode=:any OnlineFilterConfig(; grid, var_names_to_filter = ("b",),
                                                                                              velocity_names = ("u", "w"), filter_params = online_params)
end

@testset "set_offline_spectrum_filter_params" begin
    alpha, freqs = 2e-5, [1e-4, 1.4e-4]

    # In-phase and quadrature terms at each frequency, output separately, with odd sine terms
    params = set_offline_spectrum_filter_params(; alpha, freqs)
    @test params.N_coeffs == 2 && params.outputs === :separate && params.sine_parity === :odd
    @test params.a1 == params.b1 && params.a2 == params.b2
    @test params.c1 == alpha && params.c2 == alpha && params.d1 == freqs[1] && params.d2 == freqs[2]

    # Numerical integrals of the weight functions, which decay over a time 1/alpha
    τ = range(-60 / alpha, 60 / alpha, length = 2_000_001)
    integrate(f) = sum(f.(τ)) * step(τ)

    # The spectral normalisation (the default) gives the window unit energy
    @test integrate(t -> (params.a1 * exp(-alpha * abs(t)))^2) ≈ 1 rtol = 1e-4

    # The unit gain normalisation passes a signal at each frequency with its amplitude preserved
    unit_gain = set_offline_spectrum_filter_params(; alpha, freqs, normalisation = :unit_gain)
    for n in 1:2
        a, ω = getproperty(unit_gain, Symbol("a$n")), freqs[n]
        @test integrate(t -> a * exp(-alpha * abs(t)) * cos(ω * t) * cos(ω * t)) ≈ 1 rtol = 1e-4
    end
    @test set_offline_spectrum_filter_params(; alpha, freqs = [0.0], normalisation = :unit_gain).a1 ≈ alpha / 2

    @test_throws r"alpha must be positive" set_offline_spectrum_filter_params(; alpha = -1, freqs)
    @test_throws r"at least one frequency" set_offline_spectrum_filter_params(; alpha, freqs = Float64[])
    @test_throws r"normalisation must be" set_offline_spectrum_filter_params(; alpha, freqs, normalisation = :other)

    # A config with these parameters outputs the terms separately (and isn't a normalised low-pass filter),
    # so doesn't regrid
    config = @test_logs (:warn, r"setting regrid_to_mean = false") match_mode=:any OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                              var_names_to_filter = ("b",), velocity_names = ("u", "w"), filter_params = params)
    @test config.filter_params.outputs === :separate && !config.regrid_to_mean
end

@testset "set_online_spectrum_filter_params" begin
    alpha, freqs = 2e-5, [1e-4, 1.4e-4]

    # In-phase and quadrature terms at each frequency, output separately, with no sine_parity (it has no effect online)
    params = set_online_spectrum_filter_params(; alpha, freqs)
    @test params.N_coeffs == 2 && params.outputs === :separate && !haskey(params, :sine_parity)
    @test params.a1 == params.b1 && params.c1 == alpha && params.d2 == freqs[2]

    # The spectral normalisation (the default) gives the one-sided window unit energy (trapezoidal rule)
    τ = range(0, 60 / alpha, length = 1_000_001)
    energy = (params.a1 .* exp.(-alpha .* τ)).^2
    @test (sum(energy) - energy[1] / 2) * step(τ) ≈ 1 rtol = 1e-4

    # The unit gain normalisation passes a signal at each frequency through the C term with its amplitude preserved
    unit_gain = set_online_spectrum_filter_params(; alpha, freqs, normalisation = :unit_gain)
    for n in 1:2
        C_only = (; a1 = getproperty(unit_gain, Symbol("a$n")), b1 = 0.0, c1 = alpha, d1 = freqs[n], N_coeffs = 1)
        @test abs(get_frequency_response(freq = [freqs[n]], filter_params = C_only, direction = :forward)[1]) ≈ 1
    end
    @test set_online_spectrum_filter_params(; alpha, freqs = [0.0], normalisation = :unit_gain).a1 ≈ alpha

    @test_throws r"alpha must be positive" set_online_spectrum_filter_params(; alpha = -1, freqs)
    @test_throws r"normalisation must be" set_online_spectrum_filter_params(; alpha, freqs, normalisation = :other)

    # An online config with these parameters outputs the terms separately, so doesn't regrid
    grid = RectilinearGrid(size = (4, 4), x = (0, 1), z = (-1, 0), topology = (Periodic, Flat, Bounded))
    config = OnlineFilterConfig(; grid, var_names_to_filter = ("b",), velocity_names = ("u", "w"), filter_params = params)
    @test config.filter_params.outputs === :separate && !config.regrid_to_mean
end

@testset "Weight function and frequency response" begin
    zero_frequency_gain = OceananigansLagrangianFilter.Utils.zero_frequency_gain

    # The frequency response is the Fourier transform of the weight function, for each direction and sine parity
    two_terms = (; a1 = 0.7, b1 = 0.4, c1 = 0.3, d1 = 1.1, a2 = 0.2, b2 = -0.5, c2 = 0.5, d2 = 0.6, N_coeffs = 2)
    t = collect(range(-120, 120; length = 480_001)); dt = t[2] - t[1]
    freq = [0.0, 0.5, 1.1, -0.8]
    for sine_parity in (:even, :odd), direction in (:both, :forward, :backward)
        filter_params = merge(two_terms, (; sine_parity))
        G = get_weight_function(; t, tref = 0.0, filter_params, direction)
        weights = fill(dt, length(t))
        direction === :both || (weights[(length(t) + 1) ÷ 2] /= 2) # Trapezoidal rule at the jump at t = tref
        numerical = [sum(weights .* G .* exp.(-im * ω .* (0.0 .- t))) for ω in freq]
        @test isapprox(get_frequency_response(; freq, filter_params, direction), numerical; rtol = 1e-4)
    end

    # The offline response is real if the weight function is even in time, and complex otherwise
    @test get_frequency_response(; freq, filter_params = two_terms) isa Vector{Float64}
    @test get_frequency_response(; freq, filter_params = merge(two_terms, (; sine_parity = :odd))) isa Vector{ComplexF64}
    @test_throws r"direction must be" get_weight_function(; t, tref = 0.0, filter_params = two_terms, direction = :sideways)

    # The Butterworth filters are normalised, and odd sine terms don't contribute to the gain at zero frequency
    @test zero_frequency_gain(set_offline_BW2_filter_params(N = 2, freq_c = 1e-4), :both) ≈ 1
    @test zero_frequency_gain(set_offline_BW2_filter_params(N = 1, freq_c = 1e-4), :both) ≈ 1
    @test zero_frequency_gain(set_online_BW_filter_params(N = 2, freq_c = 1e-4), :forward) ≈ 1
    spectrum = set_offline_spectrum_filter_params(alpha = 1e-4, freqs = [1e-3, 2e-3])
    @test zero_frequency_gain(spectrum, :both) ≈ sum(2 * spectrum[Symbol("a$i")] * 1e-4 / (1e-8 + w^2) for (i, w) in enumerate([1e-3, 2e-3]))
end

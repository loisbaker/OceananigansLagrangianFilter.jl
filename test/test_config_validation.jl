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

    # Separate outputs can't be regridded, and don't yet work with the Eulerian filter for comparison
    separate = merge(butterworth, (outputs = :separate,))
    config = @test_logs (:warn, r"no combined fields") match_mode=:any OfflineFilterConfig(; base..., filter_params = separate)
    @test config.compute_maps && !config.regrid_to_mean
    @test_throws r"not yet supported" OfflineFilterConfig(; base..., filter_params = separate, compute_Eulerian_filter = true)

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

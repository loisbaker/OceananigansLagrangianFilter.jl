using JLD2

@testset "OfflineFilterConfig boundary relaxation validation" begin

    base = (; original_data_filename = "data/reference_sim.jld2", var_names_to_filter = ("b",),
              velocity_names = ("u", "w"), N = 1, freq_c = 1e-4)
    good_mask(x, z, p) = 1.0 # one argument per non-Flat dimension, plus relaxation_params

    # relaxation_params must contain the relaxation timescale
    @test_throws ErrorException OfflineFilterConfig(; base..., relaxation_mask_func = good_mask, relaxation_params = (; a = 1))
    @test_throws ErrorException OfflineFilterConfig(; base..., relaxation_mask_func = good_mask)

    # relaxation_mask_func must take one argument per non-Flat grid dimension (plus relaxation_params)
    bad_mask(x, p) = 1.0 # missing the z argument
    @test_throws ErrorException OfflineFilterConfig(; base..., relaxation_mask_func = bad_mask, relaxation_params = (; timescale = 1hour))

    # A fully-specified relaxation configuration constructs and stores the options
    config = OfflineFilterConfig(; base..., relaxation_mask_func = good_mask, relaxation_params = (; timescale = 1hour, a = 1))
    @test config.relaxation_mask_func === good_mask
    @test config.relaxation_params == (; timescale = 1hour, a = 1)

    # No relaxation by default
    @test isnothing(OfflineFilterConfig(; base...).relaxation_mask_func)

    # T_start + T and T_end are compared to rounding error (0.1 + 0.2 != 0.3)
    @test OfflineFilterConfig(; base..., T_start = 0.1, T = 0.2, T_end = 0.3).T_end == 0.3
    @test_throws ErrorException OfflineFilterConfig(; base..., T_start = 0.1, T = 0.2, T_end = 0.4)
end

@testset "OnlineFilterConfig boundary relaxation validation" begin

    grid = RectilinearGrid(CPU(), size = (4, 4), x = (-1, 1), z = (-1, 0),
                            topology = (Periodic, Flat, Bounded))
    base = (; grid, var_names_to_filter = ("b",), velocity_names = ("u", "w"), N = 1, freq_c = 1e-4)
    good_mask(x, z, p) = 1.0

    @test_throws ErrorException OnlineFilterConfig(; base..., relaxation_mask_func = good_mask, relaxation_params = (; a = 1))
    @test_throws ErrorException OnlineFilterConfig(; base..., relaxation_mask_func = good_mask)

    bad_mask(x, p) = 1.0
    @test_throws ErrorException OnlineFilterConfig(; base..., relaxation_mask_func = bad_mask, relaxation_params = (; timescale = 1hour))

    config = OnlineFilterConfig(; base..., relaxation_mask_func = good_mask, relaxation_params = (; timescale = 1hour, a = 1))
    @test config.relaxation_mask_func === good_mask
    @test config.relaxation_params == (; timescale = 1hour, a = 1)

    @test isnothing(OnlineFilterConfig(; base...).relaxation_mask_func)
end

@testset "Frequency scale validation" begin

    grid = RectilinearGrid(CPU(), size = (4, 4), x = (-1, 1), z = (-1, 0),
                            topology = (Periodic, Flat, Bounded))
    online = (; grid, var_names_to_filter = ("b",), velocity_names = ("u", "w"), N = 1, freq_c = 1e-4)
    offline = (; original_data_filename = "data/reference_sim.jld2", var_names_to_filter = ("b",),
                 velocity_names = ("u", "w"), N = 1, freq_c = 1e-4)

    # frequency_scale_func must take the non-Flat coordinates, the dependencies and frequency_scale_params
    @test_throws ErrorException OnlineFilterConfig(; online..., frequency_scale_func = (x, p) -> 1.0)
    @test_throws ErrorException OnlineFilterConfig(; online..., frequency_scale_func = (x, z, p) -> 1.0, frequency_scale_dependencies = (:u,))

    # Without dependencies, it must be finite and positive at every cell centre
    @test_throws ErrorException OnlineFilterConfig(; online..., frequency_scale_func = (x, z, p) -> x)

    # Dependencies need a frequency_scale_func
    @test_throws ErrorException OnlineFilterConfig(; online..., frequency_scale_dependencies = (:u,))

    # and can't be the velocities
    @test_throws ErrorException OnlineFilterConfig(; online..., frequency_scale_func = (x, z, u, p) -> 1.0, frequency_scale_dependencies = (:u,))

    # Offline, the dependencies must be read from the data
    @test_throws ErrorException OfflineFilterConfig(; offline..., frequency_scale_func = (x, z, zeta, p) -> 1.0,
                                                      frequency_scale_dependencies = (:zeta,))
    # and can't be the velocities, which are negated on the backward pass
    @test_throws ErrorException OfflineFilterConfig(; offline..., frequency_scale_func = (x, z, u, p) -> 1.0,
                                                      frequency_scale_dependencies = (:u,))
    config = OfflineFilterConfig(; offline..., frequency_scale_func = (x, z, b, p) -> 1.0 + b^2, frequency_scale_dependencies = (:b,))
    @test config.frequency_scale_dependencies == (:b,)

    # The Eulerian filter in post-processing supports a frequency scale that only depends on the coordinates, but not one
    # with dependencies (which varies in time), so that is switched off with a warning
    config = OfflineFilterConfig(; offline..., frequency_scale_func = (x, z, p) -> 1.0, compute_Eulerian_filter = true)
    @test config.compute_Eulerian_filter
    config = @test_logs (:warn, r"frequency_scale_dependencies") match_mode=:any OfflineFilterConfig(; offline...,
                            frequency_scale_func = (x, z, b, p) -> 1.0 + b^2, frequency_scale_dependencies = (:b,),
                            compute_Eulerian_filter = true)
    @test !config.compute_Eulerian_filter
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
    config = @test_logs (:warn, r"Eulerian") match_mode=:any OfflineFilterConfig(; base..., lagrangian = false)
    @test config.compute_maps && !config.regrid_to_mean && config.compute_mean_velocities

    # Likewise for un-normalised filter coefficients
    config = @test_logs (:warn, r"not normalised") match_mode=:any OfflineFilterConfig(; base..., N = nothing, freq_c = nothing,
                                                                                         filter_params = (a1 = 1.0, c1 = 1.0))
    @test config.compute_maps && !config.regrid_to_mean

    # The removed option map_to_mean gives an error explaining what to use instead
    @test_throws r"compute_maps" OfflineFilterConfig(; base..., map_to_mean = true)
end

@testset "Eulerian filter (lagrangian = false)" begin

    grid = RectilinearGrid(CPU(), size = (4, 4), x = (-1, 1), z = (-1, 0), topology = (Periodic, Flat, Bounded))
    base = (; grid, var_names_to_filter = ("b",), velocity_names = ("u", "w"), N = 1, freq_c = 1e-4)
    output_names(config) = [q.combined_name for q in OceananigansLagrangianFilter.Utils.filtered_output_quantities(config)]

    # The online config warns that the filtered variables mustn't be advected. The outputs are named as Eulerian, and
    # the mean velocities can still be computed, but not regridded
    lagrangian = OnlineFilterConfig(; base...)
    eulerian = @test_logs (:warn, r"must not be advected") match_mode=:any OnlineFilterConfig(; base..., lagrangian = false)
    @test !eulerian.regrid_to_mean && eulerian.compute_mean_velocities
    @test "b_Eulerian_filtered" in output_names(eulerian) && "u_Eulerian_filtered" in output_names(eulerian)
    @test !any(occursin("Lagrangian"), output_names(eulerian))
    @test "b_Lagrangian_filtered" in output_names(lagrangian)

    # filtered_vars_advection gives the scheme for each filtered variable, or nothing for the Eulerian filter
    @test keys(filtered_vars_advection(lagrangian)) == create_filtered_vars(lagrangian)
    @test all(scheme -> scheme isa WENO, filtered_vars_advection(lagrangian))
    @test all(scheme -> scheme isa Centered, filtered_vars_advection(lagrangian; advection = Centered()))
    @test all(isnothing, filtered_vars_advection(eulerian))

    # As a safety net, initialise_filtered_vars_from_model and create_output_fields warn if the model's advection of the
    # filtered variables doesn't match lagrangian. Only the model's advection is used, so this checks it with a stand-in
    check = OceananigansLagrangianFilter.Utils.check_filtered_vars_advection
    model_with(advection) = (; advection = merge((; momentum = WENO(), b = WENO()), advection))
    @test_logs check(model_with(filtered_vars_advection(lagrangian)), lagrangian)
    @test_logs check(model_with(filtered_vars_advection(eulerian)), eulerian)
    @test_logs (:warn, r"lagrangian = false, but the model advects") check(model_with(filtered_vars_advection(lagrangian)), eulerian)
    @test_logs (:warn, r"doesn't advect some") check(model_with(filtered_vars_advection(eulerian)), lagrangian)
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

    # The time shift (the centre of the weight function) isn't meaningful for an un-normalised filter, so there's a warning
    output_filename = tempname() * ".jld2"
    try
        jldopen(output_filename, "w") do file
            file["timeseries/t/0"] = 0.0
            file["timeseries/t/1"] = 3600.0
        end
        config = OnlineFilterConfig(; grid, var_names_to_filter = ("b",), velocity_names = ("u", "w"), filter_params = params, output_filename)
        @test_logs (:warn, r"not normalised, so the time shift") match_mode=:any compute_time_shift!(config)
        @test jldopen(file -> haskey(file, "timeseries/t_shifted/1"), output_filename)
    finally
        rm(output_filename, force = true)
    end
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

    # The responses of the separate terms add up to the whole filter's response
    for sine_parity in (:even, :odd), direction in (:both, :forward, :backward)
        filter_params = merge(two_terms, (; sine_parity))
        terms = sum(get_frequency_response(; freq, filter_params, direction, term) for term in ("C1", "S1", "C2", "S2"))
        @test terms ≈ get_frequency_response(; freq, filter_params, direction)
    end

    # A spectrum filter's cosine term is its window times the cosine, and unknown terms give an error listing the valid ones
    spectrum = set_offline_spectrum_filter_params(alpha = 0.1, freqs = [1.0, 2.0])
    @test get_weight_function(; t, tref = 0.0, filter_params = spectrum, term = "C1") ≈ spectrum.a1 .* exp.(-0.1 .* abs.(t)) .* cos.(t)
    @test_throws r"term must be one of" get_frequency_response(; freq, filter_params = spectrum, term = "C3")

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

using JLD2

@testset "Offline run test on $(arch_label(arch))" for arch in test_architectures()

    # Test the offline filter on small saved simulation data, comparing against
    # a saved reference output (generated with the same N=2 Butterworth filter).
    test_filename_stem = "data/test_offline_output_" * arch_label(arch)
    try
        filter_config = OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2", # Where the original simulation output is
                                        output_filename = test_filename_stem * ".jld2", # Where to save the filtered output
                                        var_names_to_filter = ("b",), # Which variables to filter
                                        velocity_names = ("u", "w"), # Velocities to use for Lagrangian filtering
                                        architecture = arch, # CPU(), and GPU() if CUDA finds one
                                        Δt = 20minutes, # Time step of filtering simulation
                                        T_out = 1hour, # How often to output filtered data
                                        N = 2, # Order of Butterworth filter
                                        freq_c = 1e-4 / 4, # Cut-off frequency of Butterworth filter
                                        compute_mean_velocities = true, # Whether to compute the mean velocities
                                        output_netcdf = true, # Whether to output filtered data to a netcdf file in addition to .jld2
                                        delete_intermediate_files = true, # Delete the individual output of the forward and backward passes
                                        compute_Eulerian_filter = true) # Whether to compute the Eulerian filter for comparison

        # Run the offline Lagrangian filter. Letting any error here propagate means the
        # @testset itself reports it (with a full stacktrace), rather than us swallowing
        # it into a single boolean assertion.
        run_offline_Lagrangian_filter(filter_config)

        compare_filter_output_to_reference(test_filename_stem, "data/reference_offline_output.jld2", "u")
    finally
        # Clean up generated files even if an assertion above failed.
        rm(test_filename_stem * ".jld2", force = true)
        rm(test_filename_stem * ".nc", force = true)
    end
end

@testset "Offline run test: single-exponential filter (N=1) runs without error" begin
    # The N_coeffs == 0.5 special case is exercised analytically in
    # test_initialisation.jl; here we just check that the full offline pipeline
    # (run_offline_Lagrangian_filter end-to-end) also completes for it without error.
    # There's no saved reference output for this N, so we only check the run completes
    # and produces finite, non-trivial output - not that the values are "correct".
    test_filename_stem = "data/test_offline_output_N1"
    try
        filter_config = OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                        output_filename = test_filename_stem * ".jld2",
                                        var_names_to_filter = ("b",),
                                        velocity_names = ("u", "w"),
                                        architecture = CPU(),
                                        Δt = 20minutes,
                                        T_out = 1hour,
                                        N = 1, # Single-exponential filter, instead of N=2 above
                                        freq_c = 1e-4 / 4,
                                        compute_mean_velocities = true,
                                        delete_intermediate_files = true)

        run_offline_Lagrangian_filter(filter_config)

        @test isfile(test_filename_stem * ".jld2")
        u_filtered = FieldTimeSeries(test_filename_stem * ".jld2", "u_Lagrangian_filtered")
        @test all(isfinite, u_filtered.data)
        @test any(!=(0), u_filtered.data)
    finally
        rm(test_filename_stem * ".jld2", force = true)
    end
end

@testset "Offline run test: maps output without regridding" begin
    # compute_maps = true with regrid_to_mean = false: the maps are written, but nothing is regridded
    # to the mean position, and no mean velocities are computed. Kept short (T = 3hours).
    test_filename_stem = "data/test_offline_output_maps_only"
    try
        filter_config = OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                        output_filename = test_filename_stem * ".jld2",
                                        var_names_to_filter = ("b",),
                                        velocity_names = ("u", "w"),
                                        architecture = CPU(),
                                        Δt = 20minutes,
                                        T = 3hours,
                                        T_out = 1hour,
                                        N = 2,
                                        freq_c = 1e-4 / 4,
                                        compute_maps = true,
                                        regrid_to_mean = false,
                                        compute_mean_velocities = false,
                                        delete_intermediate_files = true)

        run_offline_Lagrangian_filter(filter_config)

        output_names = jldopen(file -> keys(file["timeseries"]), test_filename_stem * ".jld2")
        @test "xi_u" in output_names && "xi_w" in output_names
        @test "b_Lagrangian_filtered" in output_names
        @test !any(endswith("_at_mean"), output_names)
        @test !("u_Lagrangian_filtered" in output_names)

        xi_u = FieldTimeSeries(test_filename_stem * ".jld2", "xi_u")
        @test all(isfinite, xi_u.data)
        @test any(!=(0), xi_u.data)
    finally
        rm(test_filename_stem * ".jld2", force = true)
    end
end

@testset "Offline run test: separate outputs and odd sine terms" begin
    # Two short runs (T = 3hours) with the same Butterworth filter and odd sine terms, output combined and separately.
    # The separate outputs of each pass, and so of the summed forward and backward passes, add up to the combined
    # outputs. This checks the backward-pass sign of the odd sine terms for both kinds of output. Only the variable
    # is filtered, to keep the runs quick (the maps and mean velocities are tested in test_initialisation.jl). The
    # Eulerian filters for comparison add up in the same way.
    stems = Dict(:combined => "data/test_offline_output_combined_odd", :separate => "data/test_offline_output_separate_odd")
    butterworth = merge(set_offline_BW2_filter_params(N = 2, freq_c = 1e-4 / 4), (sine_parity = :odd,))
    filter_params = Dict(:combined => butterworth, :separate => merge(butterworth, (outputs = :separate,)))
    try
        for (options, stem) in stems
            filter_config = OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                                output_filename = stem * ".jld2",
                                                var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                                architecture = CPU(), Δt = 20minutes, T = 3hours, T_out = 1hour,
                                                filter_params = filter_params[options],
                                                compute_maps = false, regrid_to_mean = false, compute_mean_velocities = false,
                                                compute_Eulerian_filter = true, delete_intermediate_files = true)
            run_offline_Lagrangian_filter(filter_config)
        end

        # Each pass saves its outputs as Float32, so the separate terms only add up to the combined output to Float32 precision
        output(options, name) = FieldTimeSeries(stems[options] * ".jld2", name).data
        @test isapprox(output(:separate, "b_C1_scaled") .+ output(:separate, "b_S1_scaled"),
                       output(:combined, "b_Lagrangian_filtered"); rtol = 1e-6)
        @test all(isfinite, output(:separate, "b_S1_scaled"))

        # Likewise for the Eulerian filter for comparison, which filters with each term's weight function
        @test isapprox(output(:separate, "b_Eulerian_filtered_C1_scaled") .+ output(:separate, "b_Eulerian_filtered_S1_scaled"),
                       output(:combined, "b_Eulerian_filtered"); rtol = 1e-6)
    finally
        for stem in values(stems)
            rm(stem * ".jld2", force = true)
        end
    end
end

@testset "Offline Eulerian filter in post-processing with a frequency scale" begin
    # compute_Eulerian_filter! on copies of the reference output, with a frequency scale M that only depends on the
    # coordinates. At a fixed point the weight function is then M G(M tau), so a uniform M = m is the same as a cutoff
    # m * freq_c, and a varying M is, at each point, the uniform filter with the local value of M there (evaluated where
    # each variable is saved, e.g. at the u points for u)
    freq_c = 1e-4 / 4
    config(filename; kwargs...) = OfflineFilterConfig(; original_data_filename = "data/reference_sim.jld2", output_filename = filename,
                                                      var_names_to_filter = ("b",), velocity_names = ("u", "w"), N = 2, kwargs...)
    Eulerian_filter(; kwargs...) = (filename = reference_output_without_Eulerian_filter("data/reference_offline_output.jld2");
                                    compute_Eulerian_filter!(config(filename; kwargs...)); filename)
    output(filename, name) = FieldTimeSeries(filename, name)

    m = 1.7
    scaled = Eulerian_filter(; freq_c, frequency_scale_func = (x, z, p) -> p.m, frequency_scale_params = (; m))
    faster = Eulerian_filter(; freq_c = m * freq_c)
    for name in ("b_Eulerian_filtered", "u_Eulerian_filtered", "w_Eulerian_filtered")
        @test output(scaled, name).data ≈ output(faster, name).data
    end

    M_func(x, z, p) = 1 + 0.5 * cos(2pi * x / 10kilometers) * (1 + z / 100)
    varying = Eulerian_filter(; freq_c, frequency_scale_func = M_func)
    for name in ("b_Eulerian_filtered", "u_Eulerian_filtered")
        varying_output = output(varying, name)
        @test all(isfinite, varying_output.data) # including the halo points, which are saved too
        loc = map(L -> L(), Oceananigans.Fields.location(varying_output))
        for (i, k) in ((1, 1), (6, 10)) # points at different x and z
            x, z = Oceananigans.Grids.xnode(i, 1, k, varying_output.grid, loc...), Oceananigans.Grids.znode(i, 1, k, varying_output.grid, loc...)
            uniform = Eulerian_filter(; freq_c = M_func(x, z, nothing) * freq_c)
            @test varying_output.data[i, 1, k, :] ≈ output(uniform, name).data[i, 1, k, :]
        end
    end
end

@testset "Offline Eulerian filter in post-processing reads the original data" begin
    # Offline, the Eulerian filter reads the fields from the original data, so it doesn't depend on the output interval
    # or need the unfiltered fields in the output file: with every other output time removed from the reference output
    # (as if T_out were twice the interval of the original data), and the unfiltered fields, it is the same at the
    # remaining output times
    coarse = reference_output_without_Eulerian_filter("data/reference_offline_output.jld2")
    jldopen(coarse, "a+") do file
        foreach(name -> Base.delete!(file["timeseries"], name), ("b", "u", "w"))
        iterations = sort(parse.(Int, collect(keys(file["timeseries/t"]))))
        for name in keys(file["timeseries"]), iter in iterations[2:2:end]
            Base.delete!(file, "timeseries/$name/$iter")
        end
    end
    compute_Eulerian_filter!(OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2", output_filename = coarse,
                                                 var_names_to_filter = ("b",), velocity_names = ("u", "w"), N = 2,
                                                 freq_c = 1e-4 / 4))
    for name in ("b_Eulerian_filtered", "u_Eulerian_filtered", "w_Eulerian_filtered")
        coarse_output = FieldTimeSeries(coarse, name)
        reference = FieldTimeSeries("data/reference_offline_output.jld2", name)
        @test length(coarse_output.times) < length(reference.times)
        @test coarse_output.data ≈ reference.data[:, :, :, findall(t -> t in coarse_output.times, reference.times)]
    end
end

@testset "Eulerian filter in post-processing of a constant field" begin
    # Each field is taken to be constant before the first saved time (and offline, the sums over later times start
    # from the sums over earlier times), and the weight function is integrated exactly between saved times, so a field
    # that is constant in time is unchanged by the normalised filter, including at the ends of the data
    function constant_b_copy(reference)
        filename = reference_output_without_Eulerian_filter(reference)
        jldopen(filename, "a+") do file
            for iter in filter(!=("serialized"), collect(keys(file["timeseries/b"])))
                data = file["timeseries/b/$iter"]
                Base.delete!(file, "timeseries/b/$iter")
                file["timeseries/b/$iter"] = fill(2.5, size(data))
            end
        end
        return filename
    end

    # Offline, the fields are read from the original data
    offline = reference_output_without_Eulerian_filter("data/reference_offline_output.jld2")
    compute_Eulerian_filter!(OfflineFilterConfig(original_data_filename = constant_b_copy("data/reference_sim.jld2"),
                                                 output_filename = offline, var_names_to_filter = ("b",),
                                                 velocity_names = ("u", "w"), N = 2, freq_c = 1e-4 / 4))
    @test all(FieldTimeSeries(offline, "b_Eulerian_filtered").data .≈ 2.5)

    online = constant_b_copy("data/reference_online_output.jld2")
    compute_Eulerian_filter!(OnlineFilterConfig(grid = jldopen(file -> file["serialized/grid"], online), output_filename = online,
                                                var_names_to_filter = ("b", "T"), velocity_names = ("u", "w"), N = 1,
                                                freq_c = 1e-4 / 4))
    @test all(FieldTimeSeries(online, "b_Eulerian_filtered").data .≈ 2.5)
end

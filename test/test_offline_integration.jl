@testset "Offline run test" begin

    # Test the offline filter on small saved simulation data, comparing against
    # a saved reference output (generated with the same N=2 Butterworth filter).
    test_filename_stem = "data/test_offline_output"
    try
        filter_config = OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2", # Where the original simulation output is
                                        output_filename = test_filename_stem * ".jld2", # Where to save the filtered output
                                        var_names_to_filter = ("b",), # Which variables to filter
                                        velocity_names = ("u", "w"), # Velocities to use for Lagrangian filtering
                                        architecture = CPU(), # CPU() or GPU(), if GPU() make sure you have CUDA.jl installed and imported
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

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

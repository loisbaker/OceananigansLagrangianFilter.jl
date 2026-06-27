@testset "initialise_filtered_vars_from_model" begin

    grid = RectilinearGrid(CPU(), size = (4, 4), x = (-1, 1), z = (-1, 0),
                            topology = (Periodic, Flat, Bounded))

    @testset "single-exponential filter (N_coeffs = 0.5)" begin
        filter_config = OnlineFilterConfig(grid = grid, var_names_to_filter = ("b",),
                                            velocity_names = ("u", "w"), N = 1, freq_c = 1e-4)
        filtered_vars = create_filtered_vars(filter_config)
        forcing = create_forcing(filtered_vars, filter_config)
        model = NonhydrostaticModel(grid; tracers = (filtered_vars..., :b),
                                    forcing = forcing, buoyancy = BuoyancyTracer())

        bᵢ(x, z) = 1 + x # non-uniform, so we are exercising more than just the boundary cells
        set!(model, b = bᵢ, u = 2.0) # uniform velocity so centre-interpolation is exact
        # Set w directly on the field: set!(model, w=...) would be overwritten by the
        # continuity-equation projection NonhydrostaticModel applies to enforce a
        # divergence-free velocity field (here u is uniform, so continuity forces w = 0).
        set!(model.velocities.w, 0.5)

        initialise_filtered_vars_from_model(model, filter_config)

        c1 = filter_config.filter_params.c1
        @test interior(model.tracers.b_C1) ≈ interior(model.tracers.b) ./ c1

        @test all(interior(model.tracers.xi_u_C1) .≈ (-1 / c1^2) * 2.0)
        @test all(interior(model.tracers.xi_w_C1) .≈ (-1 / c1^2) * 0.5)
    end

    @testset "multi-coefficient Butterworth filter (N_coeffs = 2)" begin
        filter_config = OnlineFilterConfig(grid = grid, var_names_to_filter = ("b",),
                                            velocity_names = ("u", "w"), N = 4, freq_c = 1e-4)
        filtered_vars = create_filtered_vars(filter_config)
        forcing = create_forcing(filtered_vars, filter_config)
        model = NonhydrostaticModel(grid; tracers = (filtered_vars..., :b),
                                    forcing = forcing, buoyancy = BuoyancyTracer())

        bᵢ(x, z) = 1 + x
        set!(model, b = bᵢ, u = 2.0) # uniform velocity so centre-interpolation is exact
        set!(model.velocities.w, 0.5) # see comment above on why w can't be set via set!(model, w=...)

        initialise_filtered_vars_from_model(model, filter_config)

        fp = filter_config.filter_params
        for i in 1:fp.N_coeffs
            ci = getproperty(fp, Symbol("c$i"))
            di = getproperty(fp, Symbol("d$i"))

            b_Ci = getproperty(model.tracers, Symbol("b_C$i"))
            b_Si = getproperty(model.tracers, Symbol("b_S$i"))
            @test interior(b_Ci) ≈ ci / (ci^2 + di^2) .* interior(model.tracers.b)
            @test interior(b_Si) ≈ di / (ci^2 + di^2) .* interior(model.tracers.b)

            xi_u_Ci = getproperty(model.tracers, Symbol("xi_u_C$i"))
            xi_u_Si = getproperty(model.tracers, Symbol("xi_u_S$i"))
            @test all(interior(xi_u_Ci) .≈ ((di^2 - ci^2) / (ci^2 + di^2)^2) * 2.0)
            @test all(interior(xi_u_Si) .≈ (-2 * ci * di / (ci^2 + di^2)^2) * 2.0)
        end
    end
end

@testset "initialise_filtered_vars_from_data" begin

    # Build the input_data NamedTuple the same way load_data does, but straight from
    # the small saved reference simulation rather than a derived "_filter_input.jld2" file.
    b_fts = FieldTimeSeries("data/reference_sim.jld2", "b")
    u_fts = FieldTimeSeries("data/reference_sim.jld2", "u")
    w_fts = FieldTimeSeries("data/reference_sim.jld2", "w")
    input_data = (var_data = (b_fts,), velocity_data = (u_fts, w_fts))

    grid = b_fts.grid
    filter_config = OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                        var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                        N = 1, freq_c = 1e-4, grid = grid)
    filtered_vars = create_filtered_vars(filter_config)
    forcing = create_forcing(filtered_vars, filter_config)
    model = NonhydrostaticModel(grid; tracers = (filtered_vars..., :b),
                                forcing = forcing, buoyancy = BuoyancyTracer())

    initialise_filtered_vars_from_data(model, input_data, filter_config)

    c1 = filter_config.filter_params.c1
    b₀ = interior(b_fts[Time(0)])
    @test interior(model.tracers.b_C1) ≈ b₀ ./ c1

    u₀_centred = interior(Field(@at (Center, Center, Center) u_fts[Time(0)]))
    w₀_centred = interior(Field(@at (Center, Center, Center) w_fts[Time(0)]))
    @test interior(model.tracers.xi_u_C1) ≈ (-1 / c1^2) .* u₀_centred
    @test interior(model.tracers.xi_w_C1) ≈ (-1 / c1^2) .* w₀_centred
end

@testset "change_sign_of_map_variables!" begin
    grid = RectilinearGrid(CPU(), size = (4, 4), x = (-1, 1), z = (-1, 0),
                            topology = (Periodic, Flat, Bounded))

    filter_config = OnlineFilterConfig(grid = grid, var_names_to_filter = ("b",),
                                        velocity_names = ("u", "w"), N = 1, freq_c = 1e-4)
    filtered_vars = create_filtered_vars(filter_config)
    forcing = create_forcing(filtered_vars, filter_config)
    model = NonhydrostaticModel(grid; tracers = (filtered_vars..., :b),
                                forcing = forcing, buoyancy = BuoyancyTracer())

    set!(model, b = 1.0, u = 2.0)
    set!(model.velocities.w, 0.5)
    initialise_filtered_vars_from_model(model, filter_config)

    xi_u_before = deepcopy(interior(model.tracers.xi_u_C1))
    xi_w_before = deepcopy(interior(model.tracers.xi_w_C1))
    b_C1_before = deepcopy(interior(model.tracers.b_C1))

    change_sign_of_map_variables!(model, filter_config)

    # The maps (xi_*) should flip sign...
    @test interior(model.tracers.xi_u_C1) ≈ -xi_u_before
    @test interior(model.tracers.xi_w_C1) ≈ -xi_w_before
    # ...but the filtered variable itself should be untouched
    @test interior(model.tracers.b_C1) ≈ b_C1_before
end

@testset "create_forcing wires in a relaxation term when boundary_relaxation = true" begin
    grid = RectilinearGrid(CPU(), size = (4, 4), x = (-1, 1), z = (-1, 0),
                            topology = (Periodic, Flat, Bounded))
    mask_func(x, z, p) = 1.0

    config_norelax = OnlineFilterConfig(grid = grid, var_names_to_filter = ("b",),
                                        velocity_names = ("u", "w"), N = 1, freq_c = 1e-4)
    config_relax = OnlineFilterConfig(grid = grid, var_names_to_filter = ("b",),
                                        velocity_names = ("u", "w"), N = 1, freq_c = 1e-4,
                                        boundary_relaxation = true, relax_timescale = 1hour,
                                        mask_func = mask_func, mask_params = (; a = 1))

    filtered_vars = create_filtered_vars(config_norelax) # same variable names regardless of relaxation
    forcing_norelax = create_forcing(filtered_vars, config_norelax)
    forcing_relax = create_forcing(filtered_vars, config_relax)

    # Without relaxation: filter forcing + original-data forcing
    @test length(forcing_norelax[:b_C1]) == 2
    @test length(forcing_norelax[:xi_u_C1]) == 2

    # With relaxation: an extra relaxation forcing term is appended
    @test length(forcing_relax[:b_C1]) == 3
    @test length(forcing_relax[:xi_u_C1]) == 3

    # Same coefficient/field-dependence structure, just with the extra term
    @test forcing_relax[:b_C1][1:2] == forcing_norelax[:b_C1]
end

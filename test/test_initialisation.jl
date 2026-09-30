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

    # Initialise from a BufferedDataReader on the small saved reference simulation. The expected
    # values come independently from FieldTimeSeries.
    b_fts = FieldTimeSeries("data/reference_sim.jld2", "b")
    u_fts = FieldTimeSeries("data/reference_sim.jld2", "u")
    w_fts = FieldTimeSeries("data/reference_sim.jld2", "w")

    grid = b_fts.grid
    filter_config = OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                        var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                        N = 1, freq_c = 1e-4, grid = grid)
    reader = create_buffered_reader(filter_config; direction = :forward)
    filtered_vars = create_filtered_vars(filter_config)
    forcing = create_forcing(filtered_vars, filter_config)
    model = NonhydrostaticModel(grid; tracers = (filtered_vars..., :b),
                                forcing = forcing, buoyancy = BuoyancyTracer())

    initialise_filtered_vars_from_data(model, reader, filter_config)

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

@testset "create_output_fields: forward and backward contributions" begin
    grid = RectilinearGrid(CPU(), size = (4, 4), x = (-1, 1), z = (-1, 0),
                            topology = (Periodic, Flat, Bounded))

    filter_config = OnlineFilterConfig(grid = grid, var_names_to_filter = ("b",),
                                        velocity_names = ("u", "w"), N = 2, freq_c = 1e-4)
    filtered_vars = create_filtered_vars(filter_config)
    forcing = create_forcing(filtered_vars, filter_config)
    model = NonhydrostaticModel(grid; tracers = (filtered_vars..., :b),
                                forcing = forcing, buoyancy = BuoyancyTracer())

    set!(model, b = (x, z) -> x + z, u = (x, z) -> 1 + x * z)
    set!(model.velocities.w, (x, z) -> x)
    initialise_filtered_vars_from_model(model, filter_config)

    forward_outputs  = create_output_fields(model, filter_config)
    backward_outputs = create_output_fields(model, filter_config; direction = :backward)
    @test keys(forward_outputs) == keys(backward_outputs)

    # The backward pass runs with the velocities negated, so its mean velocity outputs are negated
    # to give its contribution to the filtered field. The other outputs are unchanged.
    for name in keys(forward_outputs)
        forward  = interior(Field(forward_outputs[name]))
        backward = interior(Field(backward_outputs[name]))
        sign = name in ("u_Lagrangian_filtered", "w_Lagrangian_filtered") ? -1 : 1
        @test backward == sign .* forward
    end

    @test_throws ErrorException create_output_fields(model, filter_config; direction = :sideways)
end

@testset "create_output_fields: separate outputs and odd sine terms" begin
    filtered_output_names = OceananigansLagrangianFilter.Utils.filtered_output_names # internal (not exported)

    # An offline config (so that the backward pass can be tested) with a Butterworth filter, output separately
    grid = FieldTimeSeries("data/reference_sim.jld2", "b").grid
    butterworth = set_offline_BW2_filter_params(N = 2, freq_c = 1e-4)
    config(; options...) = OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                               var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                               filter_params = merge(butterworth, (; options...)), grid = grid,
                                               output_original_data = false)
    separate_even = @test_logs (:warn, r"no combined fields") match_mode=:any config(outputs = :separate)
    separate_odd  = config(outputs = :separate, sine_parity = :odd)
    combined_odd  = config(sine_parity = :odd)
    combined_even = config()

    filtered_vars = create_filtered_vars(combined_even)
    model = NonhydrostaticModel(grid; tracers = filtered_vars)
    for (n, name) in enumerate(keys(model.tracers)) # a different, non-trivial field for each tracer
        set!(getproperty(model.tracers, name), (x, z) -> n + x / 1e3 + z / 1e2)
    end
    field(output) = interior(Field(output))
    (; a1, b1, c1, d1) = butterworth

    # Separate outputs are named after their tracer with a _scaled suffix, and filtered_output_names lists them
    forward = create_output_fields(model, separate_even)
    @test sort(collect(keys(forward))) == sort(collect(filtered_output_names(separate_even)))
    @test sort(collect(keys(forward))) == sort(["b_C1_scaled", "b_S1_scaled", "xi_u_C1_scaled", "xi_u_S1_scaled",
                                                "xi_w_C1_scaled", "xi_w_S1_scaled", "u_C1_scaled", "u_S1_scaled",
                                                "w_C1_scaled", "w_S1_scaled"])
    @test field(forward["b_C1_scaled"]) == a1 .* interior(model.tracers.b_C1)
    @test field(forward["b_S1_scaled"]) == b1 .* interior(model.tracers.b_S1)
    @test field(forward["u_C1_scaled"]) ≈ -a1 .* (c1 .* interior(model.tracers.xi_u_C1) .+ d1 .* interior(model.tracers.xi_u_S1))
    @test field(forward["u_S1_scaled"]) ≈ b1 .* (d1 .* interior(model.tracers.xi_u_C1) .- c1 .* interior(model.tracers.xi_u_S1))

    # The separate outputs add up to the combined outputs
    combined = create_output_fields(model, combined_even)
    @test field(forward["b_C1_scaled"]) .+ field(forward["b_S1_scaled"]) ≈ field(combined["b_Lagrangian_filtered"])
    @test field(forward["xi_u_C1_scaled"]) .+ field(forward["xi_u_S1_scaled"]) ≈ field(combined["xi_u"])
    @test field(forward["u_C1_scaled"]) .+ field(forward["u_S1_scaled"]) ≈ field(combined["u_Lagrangian_filtered"])

    # In the backward pass, the mean velocities are negated (they change sign under time reversal), and so are
    # the sine terms if they are odd. The forward outputs don't depend on the sine parity.
    @test create_output_fields(model, separate_odd) |> keys |> collect |> sort == sort(collect(keys(forward)))
    for (options, sine_sign) in ((separate_even, 1), (separate_odd, -1))
        backward = create_output_fields(model, options; direction = :backward)
        for name in keys(forward)
            velocity_sign = startswith(name, "u_") || startswith(name, "w_") ? -1 : 1
            term_sign = occursin("_S1_", name) ? sine_sign : 1
            @test field(backward[name]) == velocity_sign * term_sign .* field(forward[name])
        end
    end

    # Combined outputs with odd sine terms: the backward contribution has its sine terms negated
    backward = create_output_fields(model, combined_odd; direction = :backward)
    @test field(backward["b_Lagrangian_filtered"]) ≈ field(forward["b_C1_scaled"]) .- field(forward["b_S1_scaled"])
    @test field(backward["u_Lagrangian_filtered"]) ≈ -(field(forward["u_C1_scaled"]) .- field(forward["u_S1_scaled"]))
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

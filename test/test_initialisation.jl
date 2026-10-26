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
    # As in the offline filter, the original variables are auxiliary fields of the model
    model = NonhydrostaticModel(grid; tracers = filtered_vars, auxiliary_fields = create_original_vars(filter_config),
                                forcing = forcing)

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

@testset "create_forcing adds a relaxation term when relaxation_mask_func is given" begin
    grid = RectilinearGrid(CPU(), size = (4, 4), x = (-1, 1), z = (-1, 0),
                            topology = (Periodic, Flat, Bounded))
    mask_func(x, z, p) = p.a + x * z # varies in space (and isn't symmetric in x and z), to check it's evaluated at the cell centre
    relax_timescale = 1hour

    for N in (1, 2) # single exponential and N_coeffs = 1
        config_norelax = OnlineFilterConfig(grid = grid, var_names_to_filter = ("b",),
                                            velocity_names = ("u", "w"), N = N, freq_c = 1e-4)
        config_relax = OnlineFilterConfig(grid = grid, var_names_to_filter = ("b",),
                                            velocity_names = ("u", "w"), N = N, freq_c = 1e-4,
                                            relaxation_mask_func = mask_func,
                                            relaxation_params = (; timescale = relax_timescale, a = 0.5))

        filtered_vars = create_filtered_vars(config_norelax) # same variable names regardless of relaxation
        model_norelax = NonhydrostaticModel(grid; tracers = (filtered_vars..., :b), forcing = create_forcing(filtered_vars, config_norelax))
        model_relax   = NonhydrostaticModel(grid; tracers = (filtered_vars..., :b), forcing = create_forcing(filtered_vars, config_relax))

        # Evaluate both forcings on the same fields: uniform velocities, so centre-interpolation is exact
        set!(model_relax, b = (x, z) -> 1 + x, u = 2.0)
        set!(model_relax.velocities.w, 0.5) # see above on why w can't be set via set!(model, w=...)
        for (n, name) in enumerate(filtered_vars)
            set!(getproperty(model_relax.tracers, name), n)
        end
        fields = merge(model_relax.velocities, model_relax.tracers, model_relax.auxiliary_fields)

        p = config_relax.filter_params
        i, k = 2, 3
        relaxation_term(name) = getproperty(model_relax.forcing, name)(i, 1, k, grid, model_relax.clock, fields) -
                                getproperty(model_norelax.forcing, name)(i, 1, k, grid, model_relax.clock, fields)
        value(name) = getproperty(fields, name)[i, 1, k]
        mask = mask_func(Oceananigans.Grids.xnode(i, grid, Center()), Oceananigans.Grids.znode(k, grid, Center()), (; a = 0.5))

        for j in 1:(p.N_coeffs == 0.5 ? 1 : p.N_coeffs)
            c = getproperty(p, Symbol("c", j))
            d = p.N_coeffs == 0.5 ? 0 : getproperty(p, Symbol("d", j))
            # Relaxation is towards the local equilibrium of each tracer (C, and S if there is one)
            targets = Dict(Symbol("b_C", j) => value(:b) * c / (c^2 + d^2),
                           Symbol("xi_u_C", j) => value(:u) * (d^2 - c^2) / (c^2 + d^2)^2,
                           Symbol("xi_w_C", j) => value(:w) * (d^2 - c^2) / (c^2 + d^2)^2)
            if p.N_coeffs != 0.5
                targets[Symbol("b_S", j)] = value(:b) * d / (c^2 + d^2)
                targets[Symbol("xi_u_S", j)] = value(:u) * (-2c * d) / (c^2 + d^2)^2
                targets[Symbol("xi_w_S", j)] = value(:w) * (-2c * d) / (c^2 + d^2)^2
            end
            for (name, target) in targets
                @test relaxation_term(name) ≈ -mask / relax_timescale * (value(name) - target)
            end
        end
    end
end

@testset "Frequency scale: forcing, initialisation and post-processing" begin
    # Each new model type takes a while to compile, so the forcing is checked on fields without a model, and the
    # checks that need a model share one, of the same type as in the create_output_fields tests above
    grid = RectilinearGrid(CPU(), size = (4, 4), x = (-1, 1), z = (-1, 0), topology = (Periodic, Flat, Bounded))
    m = 2.0
    M_func(x, z, p) = p.m * (1.5 + x) # varies in space
    M_b(x, z, b, p) = p.m * (0.5 + b) # the same, through the field b = 1 + x
    scale = (; frequency_scale_func = M_func, frequency_scale_params = (; m))
    uniform_scale = (; frequency_scale_func = (x, z, p) -> p.m, frequency_scale_params = (; m))
    relaxation = (; relaxation_mask_func = (x, z, p) -> 0.3 + x * z, relaxation_params = (; timescale = 1hour))
    config(N; freq_c = 1e-4, kwargs...) = OnlineFilterConfig(; grid, var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                                               N, freq_c, kwargs...)
    map_names(filtered_vars) = filter(name -> startswith(string(name), "xi_"), filtered_vars)

    # Fields that vary, so that each term of the forcing is different
    function set_fields!(fields, filtered_vars)
        set!(fields.b, (x, z) -> 1 + x)
        set!(fields.u, (x, z) -> 2 + z)
        set!(fields.w, 0.5)
        for (n, name) in enumerate(filtered_vars)
            set!(getproperty(fields, name), (x, z) -> n * (1 + x * z))
        end
    end
    divide_filtered_vars!(fields, filtered_vars, m) = foreach(name -> interior(getproperty(fields, name)) ./= m, filtered_vars)

    clock = Clock(time = 0.0)
    i, k = 2, 3
    x, z = Oceananigans.Grids.xnode(i, grid, Center()), Oceananigans.Grids.znode(k, grid, Center())
    M = M_func(x, z, (; m))
    for N in (1, 2) # single exponential and N_coeffs = 1
        filtered_vars = create_filtered_vars(config(N))
        fields = merge((; u = XFaceField(grid), w = ZFaceField(grid), b = CenterField(grid)),
                       NamedTuple(name => CenterField(grid) for name in filtered_vars))
        set_fields!(fields, filtered_vars)
        # The same, with the maps zero, so their forcing is just the velocity source
        fields_no_maps = merge(fields, NamedTuple(name => CenterField(grid) for name in map_names(filtered_vars)))
        forcing(F, name, f = fields) = getproperty(F, name)(i, 1, k, grid, clock, f)

        F_uniform     = create_forcing(filtered_vars, config(N))
        F_one         = create_forcing(filtered_vars, config(N; frequency_scale_func = (x, z, p) -> 1.0))
        F_scale       = create_forcing(filtered_vars, config(N; scale...))
        F_b           = create_forcing(filtered_vars, config(N; frequency_scale_func = M_b, frequency_scale_dependencies = (:b,),
                                                                frequency_scale_params = (; m)))
        F_relax       = create_forcing(filtered_vars, config(N; relaxation...))
        F_scale_relax = create_forcing(filtered_vars, config(N; scale..., relaxation...))

        for name in filtered_vars
            # M = 1 is exactly the same as no frequency scale
            @test forcing(F_one, name) == forcing(F_uniform, name)
            # A dependency on a field gives the same as the same function of the coordinates
            @test forcing(F_b, name) ≈ forcing(F_scale, name)

            relaxation_term(F_with, F_without) = forcing(F_with, name) - forcing(F_without, name)
            if name in map_names(filtered_vars)
                # Maps: the decay terms are multiplied by M, but not the velocity source...
                source = forcing(F_uniform, name, fields_no_maps)
                @test forcing(F_scale, name) ≈ source + M * (forcing(F_uniform, name) - source)
                # ...and the relaxation target (the local equilibrium) is divided by M. Without the frequency scale, the
                # relaxation term is -R/T (q - target), which gives the target
                q = getproperty(fields, name)[i, 1, k]
                R = 0.3 + x * z
                T = relaxation.relaxation_params.timescale
                target = q + relaxation_term(F_relax, F_uniform) * T / R
                @test relaxation_term(F_scale_relax, F_scale) ≈ -R / T * (q - target / M)
            else
                # Tracers: all the terms are multiplied by M, and the relaxation is unchanged
                @test forcing(F_scale, name) ≈ M * forcing(F_uniform, name)
                @test relaxation_term(F_scale_relax, F_scale) ≈ relaxation_term(F_relax, F_uniform)
            end
        end

        # A uniform frequency scale m is the same filter as a uniform filter with cutoff m * freq_c: its tracers and maps
        # are m times those of the faster filter (all the Butterworth coefficients scale with freq_c), so its forcing is too
        forcing_values(F) = Dict(name => [getproperty(F, name)(i, 1, k, grid, clock, fields) for i in 1:4, k in 1:4]
                                 for name in filtered_vars)
        forcing_scaled = forcing_values(create_forcing(filtered_vars, config(N; uniform_scale...)))
        divide_filtered_vars!(fields, filtered_vars, m)
        forcing_faster = forcing_values(create_forcing(filtered_vars, config(N; freq_c = m * 1e-4)))
        for name in filtered_vars
            @test forcing_faster[name] ≈ forcing_scaled[name] ./ m
        end
    end

    # Outputs and initialisation, for N = 2
    filtered_vars = create_filtered_vars(config(2))
    model = NonhydrostaticModel(grid; tracers = (filtered_vars..., :b))
    fields = merge(model.velocities, model.tracers)
    set_fields!(fields, filtered_vars)
    tracer(name) = interior(getproperty(model.tracers, name))
    M_field = CenterField(grid)
    set!(M_field, (x, z) -> M_func(x, z, (; m)))
    output_values(config) = Dict(name => (f = Field(output); compute!(f); copy(interior(f)))
                                 for (name, output) in create_output_fields(model, config))

    # The mean velocities are multiplied by M (to give the velocity of the mean position), and the other outputs are
    # unchanged. This uses the separate outputs, and the combined outputs are checked with the uniform frequency scale
    separate(; kwargs...) = OnlineFilterConfig(; grid, var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                                                 filter_params = merge(config(2).filter_params, (; outputs = :separate)), kwargs...)
    uniform, scaled = output_values(separate()), output_values(separate(; scale...))
    for name in keys(uniform)
        is_velocity = startswith(name, "u") || startswith(name, "w")
        @test scaled[name] ≈ (is_velocity ? interior(M_field) : 1) .* uniform[name]
    end

    # With a uniform frequency scale m, the outputs (including the mean velocities, with their factor of M) are the same
    # as those of the faster filter, with its tracers and maps 1/m times
    outputs_scaled = output_values(config(2; uniform_scale...))
    divide_filtered_vars!(fields, filtered_vars, m)
    outputs_faster = output_values(config(2; freq_c = m * 1e-4))
    for name in keys(outputs_scaled)
        @test outputs_faster[name] ≈ outputs_scaled[name]
    end

    # Initialisation: the maps' local equilibria are divided by M, and the filtered variables are unchanged...
    initialise_filtered_vars_from_model(model, config(2))
    uniform_values = Dict(name => copy(tracer(name)) for name in filtered_vars)
    initialise_filtered_vars_from_model(model, config(2; scale...))
    for name in filtered_vars
        expected = name in map_names(filtered_vars) ? uniform_values[name] ./ interior(M_field) : uniform_values[name]
        @test tracer(name) ≈ expected
    end
    # ...so with a uniform frequency scale m, they're m times those of the faster filter
    initialise_filtered_vars_from_model(model, config(2; uniform_scale...))
    scaled_values = Dict(name => copy(tracer(name)) for name in filtered_vars)
    initialise_filtered_vars_from_model(model, config(2; freq_c = m * 1e-4))
    for name in filtered_vars
        @test tracer(name) ≈ scaled_values[name] ./ m
    end

    # Dependencies that aren't fields of the model, and a frequency scale that isn't positive everywhere, are errors
    # when initialising. So are dependencies that aren't at cell centres (checked here with an extra field instead of
    # a new model, which would take a while to compile)
    @test_throws ErrorException initialise_filtered_vars_from_model(model,
        config(2; frequency_scale_func = M_b, frequency_scale_dependencies = (:vorticity,), frequency_scale_params = (; m)))
    @test_throws ErrorException initialise_filtered_vars_from_model(model,
        config(2; frequency_scale_func = (x, z, b, p) -> b - 1.5, frequency_scale_dependencies = (:b,)))
    model_with_face_field = (; velocities = model.velocities, tracers = model.tracers, auxiliary_fields = (; b_face = XFaceField(grid)))
    @test_throws ErrorException OceananigansLagrangianFilter.Utils.check_frequency_scale_dependencies(model_with_face_field,
        config(2; frequency_scale_func = M_b, frequency_scale_dependencies = (:b_face,), frequency_scale_params = (; m)))

    # The time shift and the Eulerian filter in post-processing aren't available with a frequency scale
    @test_throws ErrorException compute_time_shift!(config(1; scale...))
    # (the Eulerian filter only with frequency_scale_dependencies, which make the frequency scale vary in time)
    @test_throws ErrorException compute_Eulerian_filter!(config(1; frequency_scale_func = M_b, frequency_scale_dependencies = (:b,),
                                                                   frequency_scale_params = (; m)))
end

# Builds a small NonhydrostaticModel with online Lagrangian filtering wired in, runs it
# for a short time, and (optionally) compares the result against a saved reference.
function run_online_filter_test(filename_stem::String; N::Int, freq_c::Real)
    # Model parameters
    Nx = 10
    Nz = 10
    f = 1e-4                # Coriolis frequency [s⁻¹]
    L_front = 10kilometers  # Initial front width [m]
    aspect_ratio = 100      # L_front/H
    Ro = 0.1                # Rossby number (defines M^2)

    H = L_front / aspect_ratio  # Depth
    M² = (Ro^2 * f^2 * L_front) / H # Horizontal buoyancy gradient
    Δb = M² * L_front # Buoyancy difference across the front

    # Grid
    grid = RectilinearGrid(CPU(), size = (Nx, Nz),
                        x = (-L_front / 2, L_front / 2),
                        z = (-H, 0),
                        topology = (Periodic, Flat, Bounded))

    # Define model tracers
    tracers = (:b, :T)

    # Define model forcing
    forcing = NamedTuple()

    # Define filter configuration
    filter_config = OnlineFilterConfig(grid = grid,
                                        output_filename = filename_stem * ".jld2",
                                        var_names_to_filter = ("b", "T"),
                                        velocity_names = ("u", "w"),
                                        N = N,
                                        freq_c = freq_c)

    # Create the filtered variables - these will be tracers in the model
    filtered_vars = create_filtered_vars(filter_config)

    # Add to the existing tracers
    tracers = (filtered_vars..., tracers...)

    # Create forcing for these filtered variables
    filter_forcing = create_forcing(filtered_vars, filter_config)

    # Add these to the existing forcing
    forcing = merge(forcing, filter_forcing)

    # Define closures
    zero_filtered_var_closure = zero_closure_for_filtered_vars(filter_config)
    horizontal_closure = HorizontalScalarDiffusivity(ν = 1e-6, κ = merge((T = 1e-6, b = 1e-6), zero_filtered_var_closure))
    vertical_closure = VerticalScalarDiffusivity(ν = 1e-6, κ = merge((T = 1e-6, b = 1e-6), zero_filtered_var_closure))
    closure = (horizontal_closure, vertical_closure)

    # Define the model
    model = NonhydrostaticModel(grid;
                    coriolis = FPlane(f = f),
                    buoyancy = BuoyancyTracer(),
                    tracers = tracers,
                    forcing = forcing,
                    advection = WENO(),
                    closure = closure)

    # Initialise tracers
    # Model buoyancy and tracers
    bᵢ(x, z) = Δb * sin(2 * pi / L_front * x)
    Tᵢ(x, z) = exp(-(x / (L_front / 50)) .^ 2)
    set!(model, b = bᵢ, T = Tᵢ)

    # Set appropriate initial conditions for the filtered variables based on the actual variables
    initialise_filtered_vars_from_model(model, filter_config)

    # Define the simulation
    simulation = Simulation(model, Δt = 20minutes, stop_time = 1days)

    # Set an adaptive timestep
    conjure_time_step_wizard!(simulation, IterationInterval(20), cfl = 0.2, max_Δt = 20minutes)

    # Create filtered outputs:
    outputs = create_output_fields(model, filter_config)

    # Add in original variables if needed:
    outputs["b"] = model.tracers.b
    outputs["T"] = model.tracers.T
    outputs["u"] = model.velocities.u
    outputs["w"] = model.velocities.w

    # Output a .jld2 file:
    simulation.output_writers[:jld2fields] = JLD2Writer(
        model, outputs, filename = filter_config.output_filename, schedule = TimeInterval(1hour), overwrite_existing = true)

    # Run the simulation
    @info "Running the simulation..."

    run!(simulation)

    @info "Simulation completed in " * prettytime(simulation.run_wall_time)

    # Option to regrid to mean position
    regrid_to_mean_position!(filter_config)

    # Option to calculate Eulerian filter too
    compute_Eulerian_filter!(filter_config)

    # Option to add a shifted time coordinate
    compute_time_shift!(filter_config)

    # Option to output a final NetCDF file
    jld2_to_netcdf(filename_stem * ".jld2", filename_stem * ".nc")

    return nothing
end

@testset "Online run test" begin

    # Test the online filter on a small simulation, comparing against a saved
    # reference output (generated with the same N=1 single-exponential filter).
    filename_stem = "data/test_online_output"
    try
        # Letting any error propagate here means the @testset reports it directly
        # (with a full stacktrace), rather than us swallowing it into a boolean assertion.
        run_online_filter_test(filename_stem; N = 1, freq_c = 1e-4 / 4)

        compare_filter_output_to_reference(filename_stem, "data/reference_online_output.jld2", "w")
    finally
        rm(filename_stem * ".jld2", force = true)
        rm(filename_stem * ".nc", force = true)
    end
end

@testset "Online run test: multi-coefficient Butterworth filter (N=4) runs without error" begin
    # The N_coeffs > 0.5 branch is exercised analytically in test_initialisation.jl;
    # here we just check that the full online pipeline also completes for it without
    # error. There's no saved reference output for this N, so we only check the run
    # completes and produces finite, non-trivial output - not that the values are "correct".
    filename_stem = "data/test_online_output_N4"
    try
        run_online_filter_test(filename_stem; N = 4, freq_c = 1e-4 / 4)

        @test isfile(filename_stem * ".jld2")
        w_filtered = FieldTimeSeries(filename_stem * ".jld2", "w_Lagrangian_filtered")
        @test all(isfinite, w_filtered.data)
        @test any(!=(0), w_filtered.data)
    finally
        rm(filename_stem * ".jld2", force = true)
        rm(filename_stem * ".nc", force = true)
    end
end

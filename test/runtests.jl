# Some basic tests for OceananigansLagrangianFilter
using Test
using PythonCall
using Oceananigans.Units
using OceananigansLagrangianFilter

# TODO add more granular tests of individual functions and other configurations

@testset "OceananigansLagrangianFilter loading tests" begin
# Some trivial tests to ensure that the package and key components load correctly
    # Test that key components of OceananigansLagrangianFilter are defined
    @test isdefined(OceananigansLagrangianFilter, :OfflineFilterConfig)
    @test isdefined(OceananigansLagrangianFilter, :OnlineFilterConfig)
    @test isdefined(OceananigansLagrangianFilter, :run_offline_Lagrangian_filter)
    @test isdefined(OceananigansLagrangianFilter, :create_forcing)

    # Test that Oceananigans components are also available
    @test isdefined(OceananigansLagrangianFilter, :NonhydrostaticModel)
    
end

@testset "Offline run test" begin
    
    
    # Test the offline filter on small saved simulation data
    exception_thrown = false
    try
        # Set up the filter configuration
        filter_config = OfflineFilterConfig(original_data_filename="data/reference_sim.jld2", # Where the original simulation output is
                                        output_filename = "data/test_offline_output.jld2", # Where to save the filtered output
                                        var_names_to_filter = ("b",), # Which variables to filter
                                        velocity_names = ("u","w"), # Velocities to use for Lagrangian filtering
                                        architecture = CPU(), # CPU() or GPU(), if GPU() make sure you have CUDA.jl installed and imported
                                        Δt = 20minutes, # Time step of filtering simulation
                                        T_out = 1hour, # How often to output filtered data
                                        N = 2, # Order of Butterworth filter
                                        freq_c = 1e-4/2, # Cut-off frequency of Butterworth filter
                                        compute_mean_velocities = true, # Whether to compute the mean velocities
                                        output_netcdf = true, # Whether to output filtered data to a netcdf file in addition to .jld2
                                        delete_intermediate_files = true, # Delete the individual output of the forward and backward passes
                                        compute_Eulerian_filter = true) # Whether to compute the Eulerian filter for comparison


        # Run the offline Lagrangian filter
        run_offline_Lagrangian_filter(filter_config)
    catch e
        exception_thrown = true
        # If an error is caught, log it for debugging but allow the test framework 
        # to proceed to the final assertion.
        @warn "Code failed to execute: $e"
    end
    
    # Assert that no exception was thrown during the execution.
    @test exception_thrown == false 
    
    # We can now compare the output files to reference files stored in data/

    # Check the files have been created
    @test isfile("data/test_offline_output.jld2") 
    @test isfile("data/test_offline_output.nc")

    # Load the generated and reference data
    u_test = FieldTimeSeries("data/test_offline_output.jld2", "u")
    u_ref = FieldTimeSeries("data/reference_offline_output.jld2", "u")
    uL_test = FieldTimeSeries("data/test_offline_output.jld2", "u_Lagrangian_filtered")
    uL_ref = FieldTimeSeries("data/reference_offline_output.jld2", "u_Lagrangian_filtered")
    uLregrid_test = FieldTimeSeries("data/test_offline_output.jld2", "u_Lagrangian_filtered_at_mean")
    uLregrid_ref = FieldTimeSeries("data/reference_offline_output.jld2", "u_Lagrangian_filtered_at_mean")
    uE_test = FieldTimeSeries("data/test_offline_output.jld2", "u_Eulerian_filtered")
    uE_ref = FieldTimeSeries("data/reference_offline_output.jld2", "u_Eulerian_filtered")

    # Test vs ref original data
    @test isapprox(u_test.data, u_ref.data)

    # Test vs ref Lagrangian filtered data
    @test isapprox(uL_test.data, uL_ref.data)

    # Test vs ref Lagrangian filtered and regridded data
    @test isapprox(uLregrid_test.data, uLregrid_ref.data)

    # Test vs ref Eulerian filtered data
    @test isapprox(uE_test.data, uE_ref.data)

    # Test grids are the same
    @test u_test.grid == u_ref.grid

    # Remove test files
    rm("data/test_offline_output.jld2", force=true)
    rm("data/test_offline_output.nc", force=true)

end

@testset "Online run test" begin

    # Test the offline filter on a small simulation
    exception_thrown = false
    try
        # Model parameters
        Nx = 10
        Nz = 10
        f = 1e-4                # Coriolis frequency [s⁻¹]
        L_front = 10kilometers  # Initial front width [m]
        aspect_ratio = 100      # L_front/H
        Ro = 0.1                # Rossby number (defines M^2)


        H = L_front/aspect_ratio  # Depth
        M² = (Ro^2*f^2*L_front)/H # Horizontal buoyancy gradient
        Δb = M²*L_front # Buoyancy difference across the front


        filename_stem = "data/test_online_output";

        # Grid
        grid = RectilinearGrid(CPU(),size = (Nx, Nz), 
                            x = (-L_front/2, L_front/2),
                            z = (-H, 0),
                            topology = (Periodic, Flat, Bounded))

        # Define model tracers
        tracers = (:b,:T);

        # Define model forcing
        forcing = NamedTuple();

        # Define filter configuration
        filter_config = OnlineFilterConfig( grid = grid,
                                            output_filename = filename_stem * ".jld2",
                                            var_names_to_filter = ("b","T"), 
                                            velocity_names = ("u","w"),
                                            N = 1,
                                            freq_c = f/2)

        # Create the filtered variables - these will be tracers in the model
        filtered_vars = create_filtered_vars(filter_config)

        # Add to the existing tracers
        tracers = (filtered_vars..., tracers...)

        # Create forcing for these filtered variables
        filter_forcing = create_forcing(filtered_vars, filter_config);

        # Add these to the existing forcing
        forcing = merge(forcing, filter_forcing);

        # Define closures
        zero_filtered_var_closure = zero_closure_for_filtered_vars(filter_config)
        horizontal_closure = HorizontalScalarDiffusivity(ν=1e-6, κ=merge((T=1e-6, b= 1e-6),zero_filtered_var_closure) )
        vertical_closure = VerticalScalarDiffusivity(ν=1e-6 , κ=merge((T=1e-6, b= 1e-6),zero_filtered_var_closure) )
        closure = (horizontal_closure, vertical_closure);

        # Define the model
        model =  NonhydrostaticModel(grid;
                        coriolis = FPlane(f = f),
                        buoyancy = BuoyancyTracer(),
                        tracers = tracers,
                        forcing = forcing,
                        advection = WENO(),
                        closure = closure)

        # Initialise tracers
        # Model buoyancy and tracers
        bᵢ(x, z) = Δb*sin(2*pi/L_front * x)
        Tᵢ(x, z) = exp(-(x/(L_front/50)).^2)
        set!(model, b= bᵢ, T= Tᵢ) 

        # Set appropriate initial conditions for the filtered variables based on the actual variables
        initialise_filtered_vars_from_model(model, filter_config)

        # ### Define the simulation
        simulation = Simulation(model, Δt=20minutes, stop_time=1days)

        # ### Set an adaptive timestep
        conjure_time_step_wizard!(simulation, IterationInterval(20), cfl=0.2, max_Δt=20minutes)

        # Create filtered outputs:
        outputs = create_output_fields(model, filter_config);

        # Add in original variables if needed:
        outputs["b"] = model.tracers.b;
        outputs["T"] = model.tracers.T;
        outputs["u"] = model.velocities.u;
        outputs["w"] = model.velocities.w;

        # Output a .jld2 file:
        simulation.output_writers[:jld2fields] = JLD2Writer(
            model, outputs, filename=filter_config.output_filename, schedule=TimeInterval(1hour), overwrite_existing=true)

        # Run the simulation
        @info "Running the simulation..."

        run!(simulation)

        @info "Simulation completed in " * prettytime(simulation.run_wall_time)

        # Option to regrid to mean position
        regrid_to_mean_position!(filter_config)

        # Option to calculate Eulerian filter too
        compute_Eulerian_filter!(filter_config);

        # Option to add a shifted time coordinate
        compute_time_shift!(filter_config)

        # Option to output a final NetCDF file
        jld2_to_netcdf(filename_stem * ".jld2", filename_stem * ".nc")
    catch e
        exception_thrown = true
        # If an error is caught, log it for debugging but allow the test framework 
        # to proceed to the final assertion.
        @warn "Code failed to execute: $e"
    end

    # Assert that no exception was thrown during the execution.
    @test exception_thrown == false 
    
    # We can now compare the output files to reference files stored in data/

    # Check the files have been created
    @test isfile("data/test_online_output.jld2") 
    @test isfile("data/test_online_output.nc")

    # Load the generated and reference data
    w_test = FieldTimeSeries("data/test_online_output.jld2", "w")
    w_ref = FieldTimeSeries("data/reference_online_output.jld2", "w")
    wL_test = FieldTimeSeries("data/test_online_output.jld2", "w_Lagrangian_filtered")
    wL_ref = FieldTimeSeries("data/reference_online_output.jld2", "w_Lagrangian_filtered")
    wLregrid_test = FieldTimeSeries("data/test_online_output.jld2", "w_Lagrangian_filtered_at_mean")
    wLregrid_ref = FieldTimeSeries("data/reference_online_output.jld2", "w_Lagrangian_filtered_at_mean")
    wE_test = FieldTimeSeries("data/test_online_output.jld2", "w_Eulerian_filtered")
    wE_ref = FieldTimeSeries("data/reference_online_output.jld2", "w_Eulerian_filtered")

    # Test vs ref original data
    @test isapprox(w_test.data, w_ref.data)

    # Test vs ref Lagrangian filtered data
    @test isapprox(wL_test.data, wL_ref.data)

    # Test vs ref Lagrangian filtered and regridded data
    @test isapprox(wLregrid_test.data, wLregrid_ref.data)

    # Test vs ref Eulerian filtered data
    @test isapprox(wE_test.data, wE_ref.data)

    # Test grids are the same
    @test w_test.grid == w_ref.grid


    # Remove test files
    rm("data/test_online_output.jld2", force=true)
    rm("data/test_online_output.nc", force=true)
end
using Oceananigans.TimeSteppers: reset!
using Printf

# Run the filtering operation

function run_offline_Lagrangian_filter(config)

    # Copy and manipulate data on disk to have correct order and time shift
    create_input_data_on_disk(config; direction = "forward") 

    # Load in saved data from simulation
    saved_velocities, saved_original_vars = load_data(config)
    println("Loaded data from $(config.original_data_filename)")

    # Create the original variables - these will be auxiliary fields in the model
    original_vars = create_original_vars(config)
    println("Created original variables: ", original_vars)

    # Create the filtered variables - these will be tracers in the model
    filtered_vars = create_filtered_vars(config)
    println("Created filtered variables: ", filtered_vars)

    # Create forcing for these filtered variables
    forcing = create_forcing(filtered_vars, config)
    println("Created forcing for filtered variables ")

    # Define model 
    model = LagrangianFilter(;config.grid, tracers = filtered_vars, auxiliary_fields = original_vars, forcing = forcing, advection=config.advection)
    println("Created model")

    # We can set initial values to improve the spinup, use the limit freq_c -> \infty
    # The map variables get automatically initialised to zero
    initialise_filtered_vars(model, saved_original_vars, config)    
    println("Initialised filtered variables")

    # Define our outputs # 
    filtered_outputs = create_output_fields(model, config)
    println("Defined outputs")

    # Define the filtering simulation 
    simulation = Simulation(model, Δt = config.Δt, stop_time = config.T) 
    println("Defined simulation")

    # Tell the simulation to use the saved data
    simulation.callbacks[:update_input_data] = Callback(update_input_data!, parameters = (velocities = saved_velocities, original_vars = saved_original_vars))

    # Add a progress monitor
    function progress(sim)
        @info @sprintf("Simulation time: %s\n", 
                    prettytime(sim.model.clock.time))             
        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress,TimeInterval(config.T_out))

    #Write outputs
    simulation.output_writers[:vars] = JLD2Writer(model, filtered_outputs,
                                                            filename = config.forward_output_filename,
                                                            schedule = TimeInterval(config.T_out),
                                                            overwrite_existing = true)

    # Run forward simulation                                                        
    run!(simulation)


    # Now, run it backwards. Switch the data direction on disk
    create_input_data_on_disk(config; direction = "backward")

    # Reload the saved data
    saved_velocities, saved_original_vars = load_data(config)

    # Tracers are initialised with their existing values 
  
    # Reset time
    reset!(model.clock)

    # Write outputs
    simulation.output_writers[:vars] = JLD2Writer(model, filtered_outputs,
                                                            filename = config.backward_output_filename,
                                                            schedule = TimeInterval(config.T_out),
                                                            overwrite_existing = true)

    # And run the backward simulation.
    run!(simulation)

    # Now sum the forward and backward components
    sum_forward_backward_contributions!(config)

    # Clean up temporary files
    if config.delete_intermediate_files
        rm(config.forward_output_filename)
        rm(config.backward_output_filename)
        rm(config.original_data_filename[1:end-5] * "_filter_input.jld2")
    end

    # Regrid if necessary
    if config.map_to_mean
        regrid_to_mean_position!(config)
    end

    # Output netcdf if necessary
    if config.output_netcdf
        jld2_to_netcdf(config.combined_output_filename, config.combined_output_filename[1:end-5] * ".nc")
    end

end

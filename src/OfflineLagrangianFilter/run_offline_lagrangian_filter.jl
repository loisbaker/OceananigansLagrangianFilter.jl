using Oceananigans.TimeSteppers: reset!
using Printf

# Run the filtering operation. This is the main function that a user will call to perform offline filtering.
# The utility functions used here can also be called individually if desired to create a custom filtering workflow.

"""
    run_offline_Lagrangian_filter(config::OfflineFilterConfig)

Runs an offline Lagrangian filter on an Oceananigans `FieldTimeSeries` dataset as configured by `config`.

This function performs a series of steps to filter the data:
1.  **Run forward simulation**: A `LagrangianFilter` model is created and run forward in time to compute the first half of the filter contributions. Input data is read directly from the original source file via a `BufferedDataReader` — no intermediate file is written.
2.  **Run backward simulation**: A second `BufferedDataReader` is created for the backward pass and the simulation is run again. Velocity reversal is handled on-the-fly during interpolation.
3.  **Combine results**: The forward and backward simulation outputs are summed to produce the final filtered data.
4.  **Post-processing**: Optional post-processing steps are performed, including regridding the data to the mean position, computing a comparative Eulerian filter, and converting the output file to NetCDF.

Arguments
=========

- `config`: An instance of `OfflineFilterConfig` that specifies all parameters and file paths for the filtering process.
"""
function run_offline_Lagrangian_filter(config)

    # Create the original variables - these will be auxiliary fields in the model
    original_vars = create_original_vars(config)
    @info "Created original variables: $(keys(original_vars))"

    # Create the filtered variables - these will be tracers in the model
    filtered_vars = create_filtered_vars(config)
    @info "Created filtered variables: $filtered_vars"

    # Create forcing for these filtered variables
    forcing = create_forcing(filtered_vars, config)
    @info "Created forcing for filtered variables"

    # Define model
    model = LagrangianFilter(config.grid; tracers = filtered_vars, auxiliary_fields = original_vars, forcing = forcing, advection=config.advection)
    @info "Created model"

    # ── Forward pass ──────────────────────────────────────────────────────────

    # Open a buffered reader for the forward direction — reads directly from the
    # source file with a two-frame GPU buffer; no intermediate file is written.
    forward_reader = create_buffered_reader(config; direction = :forward)
    @info "Loaded forward input data from $(config.original_data_filename)"

    # Initialise filtered variables from the data at t=0
    initialise_filtered_vars_from_data(model, forward_reader, config)
    @info "Initialised filtered variables"

    # Define our outputs
    filtered_outputs = create_output_fields(model, config)
    @info "Defined outputs"

    # Define the filtering simulation
    simulation = Simulation(model, Δt = config.Δt, stop_time = config.T)
    @info "Defined simulation"

    # Tell the simulation to use the buffered reader.
    simulation.callbacks[:update_input_data] = Callback(update_input_data!, callsite = UpdateStateCallsite(), parameters = forward_reader)

    # Add a progress monitor
    function progress(sim)
        @info @sprintf("Simulation time: %s\n",
                    prettytime(sim.model.clock.time))
        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, TimeInterval(config.T/10))

    # Write outputs
    simulation.output_writers[:vars] = JLD2Writer(model, filtered_outputs,
                                                            filename = config.forward_output_filename,
                                                            schedule = TimeInterval(config.T_out),
                                                            overwrite_existing = true)

    # Run forward simulation
    run!(simulation)

    # ── Backward pass ─────────────────────────────────────────────────────────

    # Create a new buffered reader for the backward direction.
    # Velocity negation for backward advection is handled inside interpolate_to_model!.
    backward_reader = create_buffered_reader(config; direction = :backward)
    @info "Loaded backward input data from $(config.original_data_filename)"

    # The filtered variables are already well initialised from the forward run,
    # but the map variables need their sign reversed.
    change_sign_of_map_variables!(model, config)

    # Reset simulation time to zero
    reset!(model.clock)

    # Swap in the backward reader
    simulation.callbacks[:update_input_data] = Callback(update_input_data!, callsite = UpdateStateCallsite(), parameters = backward_reader)

    # Write outputs
    simulation.output_writers[:vars] = JLD2Writer(model, filtered_outputs,
                                                            filename = config.backward_output_filename,
                                                            schedule = TimeInterval(config.T_out),
                                                            overwrite_existing = true)

    # Run the backward simulation
    run!(simulation)

    # ── Post-processing ───────────────────────────────────────────────────────

    # Sum the forward and backward components
    sum_forward_backward_contributions!(config)

    # Option to delete the forward and backward output files
    if config.delete_intermediate_files
        rm(config.forward_output_filename)
        rm(config.backward_output_filename)
    end

    # Option to regrid to mean position
    if config.map_to_mean
        regrid_to_mean_position!(config)
    end

    # Option to calculate Eulerian filter too
    if config.compute_Eulerian_filter
        compute_Eulerian_filter!(config)
    end

    # Option to output final netcdf
    if config.output_netcdf
        jld2_to_netcdf(config.output_filename, config.output_filename[1:end-5] * ".nc")
    end

end

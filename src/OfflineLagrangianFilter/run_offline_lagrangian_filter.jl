using Oceananigans.TimeSteppers: reset!
using Printf

# Run the filtering operation. This is the main function that a user will call to perform offline filtering.
# The steps of the filter itself are in run_forward_and_backward_passes: the utility functions used there can also be
# called individually if desired to create a custom filtering workflow. Running the filter in batches of terms
# (batch_size) is in filter_term_batches.jl.

"""
    run_offline_Lagrangian_filter(config::OfflineFilterConfig; batch_size::Union{Int, Nothing} = nothing)

Runs an offline Lagrangian filter on an Oceananigans `FieldTimeSeries` dataset as configured by `config`.

This function performs a series of steps to filter the data:
1.  **Run forward simulation**: A `LagrangianFilter` model is created and run forward in time to compute the first half of the filter contributions. Input data is read directly from the original source file via a `BufferedDataReader` — no intermediate file is written.
2.  **Run backward simulation**: A second `BufferedDataReader` is created for the backward pass and the simulation is run again. Velocity reversal is handled on-the-fly during interpolation.
3.  **Combine results**: The forward and backward simulation outputs are summed to produce the final filtered data.
4.  **Post-processing**: Optional post-processing steps are performed, including regridding the data to the mean position, computing a comparative Eulerian filter, and converting the output file to NetCDF.

The filter can also be run in batches of `batch_size` terms (each term of the filter, e.g. each frequency of a
spectrum filter, has a cosine and a sine tracer for each filtered variable and map), one after the other. Each batch
runs steps 1-3 with a smaller model, writing to its own file (`<output_filename>_batch<n>.jld2`), which is merged into
`output_filename` as soon as the batch is done and then deleted (unless `delete_intermediate_files = false`).
Outputs of separate terms (`outputs = :separate`) are numbered as in the whole filter, and combined outputs are the
sums of those of the batches. The output is then the same as without batches (combined outputs to the Float32
precision they are saved with), and the post-processing in step 4 is done once on the merged output. Batches of the
same size have the same model types, so only the first is compiled (and the last, if it is smaller). Batches help when
a filter has many terms, as the compilation time grows faster than the number of terms, and a smaller model needs less
memory.

Arguments
=========

- `config`: An instance of `OfflineFilterConfig` that specifies all parameters and file paths for the filtering process.
- `batch_size`: The number of filter terms in each batch, or `nothing` (default) to run all the terms together.
"""
function run_offline_Lagrangian_filter(config; batch_size::Union{Int, Nothing} = nothing)

    N_terms = config.filter_params.N_coeffs
    if isnothing(batch_size) || N_terms == 0.5 || batch_size >= N_terms # N_terms = 0.5 is a single exponential
        run_forward_and_backward_passes(config)
    else
        run_batches(config, batch_size)
    end

    # ── Post-processing ───────────────────────────────────────────────────────

    # Option to regrid to mean position
    if config.regrid_to_mean
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

# Run the forward and backward passes of the filter in `config`, and sum them into `config.output_filename`
function run_forward_and_backward_passes(config)

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
    model = LagrangianFilter(config.grid; tracers = filtered_vars, auxiliary_fields = original_vars, forcing = forcing, advection = config.lagrangian ? config.advection : nothing)
    @info "Created model"

    # ── Forward pass ──────────────────────────────────────────────────────────

    # Open a buffered reader for the forward direction — reads directly from the
    # source file with a two-frame GPU buffer; no intermediate file is written.
    forward_reader = create_buffered_reader(config; direction = :forward)
    @info "Loaded forward input data from $(config.original_data_filename)"

    # Initialise filtered variables from the data at t=0
    initialise_filtered_vars_from_data(model, forward_reader, config)
    @info "Initialised filtered variables"

    # Define our outputs. Oceananigans names interpolation operators with a global counter, so building the same
    # outputs again for the backward pass would give new types, and the output writer would be compiled again.
    # Restoring the counter before the backward pass gives the same types (if a future Oceananigans has no such
    # counter, the backward pass is just compiled again).
    identity_counter = output_identity_counter()
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
                                                            overwrite_files = true)

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

    # Write outputs: the backward pass's contributions to the filtered fields
    restore_output_identity_counter!(identity_counter)
    backward_outputs = create_output_fields(model, config; direction = :backward)
    simulation.output_writers[:vars] = JLD2Writer(model, backward_outputs,
                                                            filename = config.backward_output_filename,
                                                            schedule = TimeInterval(config.T_out),
                                                            overwrite_files = true)

    # Run the backward simulation
    run!(simulation)

    # Sum the forward and backward components
    sum_forward_backward_contributions!(config)

    # Option to delete the forward and backward output files
    if config.delete_intermediate_files
        rm(config.forward_output_filename)
        rm(config.backward_output_filename)
    end

    return nothing
end

# The global counter that Oceananigans uses to name interpolation operators, or nothing if there is none
output_identity_counter() = isdefined(Oceananigans.Operators, :identity_counter) ? Oceananigans.Operators.identity_counter : nothing

restore_output_identity_counter!(::Nothing) = nothing
restore_output_identity_counter!(counter) = setglobal!(Oceananigans.Operators, :identity_counter, counter)

using OceananigansLagrangianFilter
using Printf
using NCDatasets

# User defined options

fields_filename = joinpath(@__DIR__, "SW_vort.jld2")
T_start = 0
T_end = 4

arch = GPU()

# Set the output period
T_out = 0.05

# Set filter order and cut-off frequency
# Amplitude of frequency response of filter will be squared Butterworth order 2^N
N = 2 
freq_c = 2

# Define tracers to filter
filter_tracer_names = ("ω",)

# Define velocities to use for filtering
velocity_names = ("u","v")

# Set filtering parameters (this is for Butterworth-type, could define others here)
filter_params = set_BW_filter_params(N=N,freq_c=freq_c)

# Set the time step for the simulation
Δt = 1e-3

# Decide whether to solve for and output maps to generalised Lagrangian mean
map_to_mean = true

# Name output files
forward_output_filename = joinpath(@__DIR__, "forward_LF.jld2")
backward_output_filename = joinpath(@__DIR__, "backward_LF.jld2")
combined_output_filename = joinpath(@__DIR__, "combined_LF.jld2")

# Manipulate data on disk to have correct order 
T = set_data_on_disk!(fields_filename, direction="forward", T_start = T_start, T_end = T_end)

# Load in saved data from simulation
saved_velocities, saved_tracers, grid = load_data(fields_filename, filter_tracer_names, velocity_names, architecture=arch, backend=InMemory(4))
println("Data loaded")
# Create all the tracers we'll need to solve for
tracers = create_tracers(filter_tracer_names, velocity_names, filter_params, map_to_mean=map_to_mean)
println("Tracers created")
# Create forcing for these tracers
forcing = create_forcing(tracers, saved_tracers, filter_tracer_names, velocity_names, filter_params)
println("Forcing created")
# Define model 
model = LagrangianFilter(;grid, tracers = tracers, forcing = forcing)
println("Model created")
# Define our outputs
filtered_outputs = create_output_fields(model, filter_tracer_names, velocity_names, filter_params)
println("Outputs created")
# Define the filtering simulation 
simulation = Simulation(model, Δt = Δt, stop_time = T) 
println("Simulation created")
# Tell the simulation to use the saved velocities
simulation.callbacks[:update_velocities] = Callback(update_velocities!, parameters = saved_velocities)

# Add a progress monitor
function progress(sim)
    @info @sprintf("Simulation time: %s, max(|u|), min(|u|): %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, model.velocities.u), minimum(abs, model.velocities.u))             
     return nothing
 end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

# Write outputs
simulation.output_writers[:fieldsjld2] = JLD2Writer(model, filtered_outputs,
                                                        filename = forward_output_filename,
                                                        schedule = TimeInterval(T_out),
                                                        overwrite_existing = true)

                                                        # Write outputs
simulation.output_writers[:fieldsnc] = NetCDFWriter(model, filtered_outputs,
                                                        filename = joinpath(@__DIR__, "forward_LF.nc"),
                                                        schedule = TimeInterval(T_out),
                                                        overwrite_existing = true)

simulation.output_writers[:fieldsncall] = NetCDFWriter(model, (ωC1=model.tracers.ωC1, ωC2=model.tracers.ωC2, ωS1=model.tracers.ωS1, ωS2=model.tracers.ωS2,u=model.velocities.u,v=model.velocities.v),
                                                        filename = joinpath(@__DIR__, "all_tracers_forward.nc"),
                                                        schedule = TimeInterval(T_out),
                                                        overwrite_existing = true)
# And run the forward simulation
run!(simulation)
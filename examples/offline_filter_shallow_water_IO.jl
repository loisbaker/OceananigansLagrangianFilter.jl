# # Shallow water intertial oscillations with offline Lagrangian filtering

# This example demonstrates how to perform offline filtering on a shallow water simulation to
# remove the effect of inertial oscillations on a tracer field.

# We could also filter a shallow water simulation online, but would have to use the 
# VectorInvariantFormulation() in order to have direct access to the model velocities. 
# This example uses the ConservativeFormulation() instead, and filtering is performed 
# offline using the saved velocities after they have been calculated from uh and vh.

# ## Run the simulation 
# ### Install dependencies

using Oceananigans
using Printf
using NCDatasets

filename_stem = "SW_IO_with_tracer";

# ### Define the grid

grid = RectilinearGrid(CPU(), size = (50, 50),
                       x = (0, 2*pi),
                       y = (0, 2*pi),
                       topology = (Periodic, Periodic, Flat))

# ### Set parameters
# Building a `ShallowWaterModel`. We non-dimensionalise as in Kafiabad & Vanneste 2023.
Fr = 0.1 # Froude number
Ro = 1 # fRossby number

gravitational_acceleration = 1/Fr^2
coriolis = FPlane(f=1/Ro)

# ### Define the model
model = ShallowWaterModel(; grid, coriolis, gravitational_acceleration,
                            timestepper = :RungeKutta3,
                            tracers= (:T,),
                            momentum_advection = WENO())


# ### Initial conditions

# Velocity and height initial conditions - uniform velocity perturbation, initial height is 1 (unperturbed)
displacement = 2*pi/10
u_i = displacement/Ro   
h_i = 1
uh_i = u_i*h_i;

# Initialise a tracer as a blob in the middle of the domain
width = 2*pi/15
T_i(x, y) = exp(-(((x - pi)^2 + (y - pi)^2)/width).^2)

set!(model, uh = uh_i, h= h_i, T = T_i )

uh, vh, h = model.solution

u = Field(uh / h)
v = Field(vh / h)
T = model.tracers.T

# ### Simulation
simulation = Simulation(model, Δt = 1e-2, stop_time = 20)

function progress(sim)
    model = sim.model
    uh, vh, h = model.solution
    @info @sprintf("Simulation time: %s, max(|uh|, |vh|, |h|): %.2e, %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, uh), maximum(abs, vh), 
                   maximum(abs, h))

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# ### Set up the outputs
# Save velocities and tracer for Lagrangian filtering
simulation.output_writers[:fields_jld2] = JLD2Writer(model, (; u,v,T),
                                                        filename = filename_stem * ".jld2",
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)

# ### And finally run the simulation.
run!(simulation)

# ## Perform Lagrangian filtering
# Now we set up and run the offline Lagrangian filter on the output of the above simulation.
# This could be performed in a different script (with appropriate import of Oceananigans.Units and CUDA if needed)

using OceananigansLagrangianFilter

# ### Set up the filter configuration
filter_config = OfflineFilterConfig(original_data_filename="SW_IO_with_tracer.jld2", # Where the original simulation output is
                                    output_filename = "SW_IO_with_tracer_filtered.jld2", # Where to save the filtered output
                                    var_names_to_filter = ("T",), # Which variables to filter
                                    velocity_names = ("u","v"), # Velocities to use for Lagrangian filtering
                                    architecture = CPU(), # CPU() or GPU(), if GPU() make sure you have CUDA.jl installed and imported
                                    Δt = 1e-2, # Time step of filtering simulation
                                    T_out=0.1, # How often to output filtered data
                                    N=2, # Order of Butterworth filter
                                    freq_c = 0.5, # Cut-off frequency of Butterworth filter
                                    compute_mean_velocities= true, # Whether to compute mean velocities
                                    output_netcdf = true, # Whether to output filtered data to a netcdf file in addition to .jld2
                                    delete_intermediate_files = true, # Delete the individual output of the forward and backward passes
                                    compute_Eulerian_filter = true) # Whether to compute the Eulerian filter for comparison


# ### Run the offline Lagrangian filter
run_offline_Lagrangian_filter(filter_config)

# ### Visualisation

using CairoMakie 

# Now we animate the results. 
timeseries1 = FieldTimeSeries(filter_config.output_filename, "T")
timeseries2 = FieldTimeSeries(filter_config.output_filename, "T_Eulerian_filtered")
timeseries3 = FieldTimeSeries(filter_config.output_filename, "T_Lagrangian_filtered")
timeseries4 = FieldTimeSeries(filter_config.output_filename, "T_Lagrangian_filtered_at_mean")

times = timeseries1.times

set_theme!(Theme(fontsize = 20))
fig = Figure(size = (1300, 500))

axis_kwargs = (xlabel = "x",
               ylabel = "y",
               limits = ((0, 2*pi), (0, 2*pi)),
               aspect = AxisAspect(1))

ax1 = Axis(fig[2, 1]; title = "Raw", axis_kwargs...)
ax2 = Axis(fig[2, 2]; title = "Eulerian filtered", axis_kwargs...)
ax3 = Axis(fig[2, 3]; title = "Lagrangian filtered", axis_kwargs...)
ax4 = Axis(fig[2, 4]; title = "Lagrangian filtered at mean", axis_kwargs...)


n = Observable(1)
Observable(1)

var1 = @lift timeseries1[$n]
var2 = @lift timeseries2[$n]
var3 = @lift timeseries3[$n]
var4 = @lift timeseries4[$n]

heatmap!(ax1, var1; colormap = :Spectral, colorrange = (0, 1))
heatmap!(ax2, var2; colormap = :Spectral, colorrange = (0, 1))
heatmap!(ax3, var3; colormap = :Spectral, colorrange = (0, 1))
heatmap!(ax4, var4; colormap = :Spectral, colorrange = (0, 1))


title = @lift "Tracer T at time = " * string(round(times[$n], digits=2))
Label(fig[1, 1:4], title, fontsize=24, tellwidth=false)

fig

frames = 1:length(times)

@info "Making an animation"

CairoMakie.record(fig, "IO_filtered_tracer_movie_offline.mp4", frames, framerate=24) do i
    n[] = i
end

# ![](IO_filtered_tracer_movie_offline.mp4)

# We see that the Eulerian filter smudges the tracer field as it is advected by the 
# inertial oscillations, while the Lagrangian filter is able to effectively remove 
# the oscillations while preserving the tracer structures.

# We remove these files to keep things tidy, keep them for analysis if desired
rm(filename_stem * ".jld2")
rm(filter_config.output_filename)
rm(filter_config.output_filename[1:end-5] * ".nc")
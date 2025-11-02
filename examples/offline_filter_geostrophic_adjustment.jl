# # Geostrophic adjustment with offline Lagrangian filtering

# We set up a geostrophic adjustment problem similar to Blumen (2000), JPO
# in a domain that is horizontally periodic. 

# An initially unbalanced two-dimensional front oscillates with the inertial 
# frequency around a state of geostrophic balance, and we illustrate that we 
# can remove the oscillations to find the mean state.
# Credit to Tom Cummings for work on this example. 

# In this example, the filtering is performed offline after the simulation.

# ## Run the simulation
# ### Install dependencies

using Oceananigans 
using Oceananigans.Units
using NCDatasets
using Printf


# ### Model parameters
Nx = 400
Nz = 80
f = 1e-4                # Coriolis frequency [s⁻¹]
L_front = 10kilometers  # Initial front width [m]
aspect_ratio = 100      # L_front/H
Ro = 0.1                # Rossby number (defines M^2)

H = L_front/aspect_ratio  # Depth
M² = (Ro^2*f^2*L_front)/H # Horizontal buoyancy gradient
Δb = M²*L_front # Buoyancy difference across the front
κh = 1e-6 # Horizontal diffusivity
κv = 1e-6 # Vertical diffusivity

filename_stem = "geostrophic_adjustment";

# ### Define the grid

grid = RectilinearGrid(CPU(),size = (Nx, Nz), 
                       x = (-L_front/2, L_front/2),
                       z = (-H, 0),
                       topology = (Periodic, Flat, Bounded))

# ### Closures                     
horizontal_closure = HorizontalScalarDiffusivity(ν=κh, κ=κh )
vertical_closure = VerticalScalarDiffusivity(ν=κv , κ=κv )
closure = (horizontal_closure, vertical_closure);

# ### Tracers
tracers = (:b,:T)

# ### Define the model
model =  NonhydrostaticModel(; grid,
                coriolis = FPlane(f = f),
                buoyancy = BuoyancyTracer(),
                tracers = tracers,
                advection = WENO(),
                closure = closure)


# ### Initialise the buoyancy and tracer 
bᵢ(x, z) = Δb*sin(2*pi/L_front * x)
Tᵢ(x, z) = exp(-(x/(L_front/50)).^2)
set!(model, b= bᵢ, T= Tᵢ) 

# ### Define the simulation
simulation = Simulation(model, Δt=20minutes, stop_time=3days)

# Set an adaptive timestep
conjure_time_step_wizard!(simulation, IterationInterval(20), cfl=0.2, max_Δt=20minutes)

# Add a progress callback

wall_clock = Ref(time_ns())

function print_progress(sim)
    u, v, w = model.velocities
    progress = 100 * (time(sim) / sim.stop_time)
    elapsed = (time_ns() - wall_clock[]) / 1e9

    @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u): (%6.3e, %6.3e, %6.3e) m/s, next Δt: %s\n",
            progress, iteration(sim), prettytime(sim), prettytime(elapsed),
            maximum(abs, u), maximum(abs, v), maximum(abs, w), prettytime(sim.Δt))

    wall_clock[] = time_ns()

    return nothing
end

add_callback!(simulation, print_progress, IterationInterval(50))

# ### Set up the output 
u, v, w = model.velocities
b = model.tracers.b
T = model.tracers.T

# Output a jld2 file for Lagrangian filtering
simulation.output_writers[:jld2fields] = JLD2Writer(
    model, (; b, u, v, w, T), filename = filename_stem * ".jld2", schedule=TimeInterval(1hour), overwrite_existing=true)

# ### Run simulation
@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)

# ## Perform Lagrangian filtering
# Now we set up and run the offline Lagrangian filter on the output of the above simulation.
# This could be performed in a different script (with appropriate import of Oceananigans.Units and CUDA if needed)

using OceananigansLagrangianFilter

# ### Set up the filter configuration
filter_config = OfflineFilterConfig(original_data_filename="geostrophic_adjustment.jld2", # Where the original simulation output is
                                    output_filename = "geostrophic_adjustment_offline_filtered.jld2", # Where to save the filtered output
                                    var_names_to_filter = ("T", "b"), # Which variables to filter
                                    velocity_names = ("u","w","v"), # Velocities to use for remapping
                                    architecture = CPU(), # CPU() or GPU(), if GPU() make sure you have CUDA.jl installed and imported
                                    Δt = 20minutes, # Time step of filtering simulation
                                    T_out = 1hour, # How often to output filtered data
                                    N = 2, # Order of Butterworth filter
                                    freq_c = 1e-4/2, # Cut-off frequency of Butterworth filter
                                    compute_mean_velocities = true, # Whether to compute the mean velocities
                                    output_netcdf = true, # Whether to output filtered data to a netcdf file in addition to .jld2
                                    delete_intermediate_files = true, # Delete the individual output of the forward and backward passes
                                    compute_Eulerian_filter = true) # Whether to compute the Eulerian filter for comparison


# ### Run the offline Lagrangian filter
run_offline_Lagrangian_filter(filter_config)

# ### Visualisation

using CairoMakie 

# Now we animate the results. First, the buoyancy:
timeseries1 = FieldTimeSeries(filter_config.output_filename, "b")
timeseries2 = FieldTimeSeries(filter_config.output_filename, "b_Eulerian_filtered")
timeseries3 = FieldTimeSeries(filter_config.output_filename, "b_Lagrangian_filtered")
timeseries4 = FieldTimeSeries(filter_config.output_filename, "b_Lagrangian_filtered_at_mean")

times = timeseries1.times

set_theme!(Theme(fontsize = 20))
fig = Figure(size = (1300, 500))

axis_kwargs = (xlabel = "x",
               ylabel = "z",
               limits = ((-5000, 5000), (-100, 0)),
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

heatmap!(ax1, var1; colormap = :balance, colorrange = (-1e-4, 1e-4))
heatmap!(ax2, var2; colormap = :balance, colorrange = (-1e-4, 1e-4))
heatmap!(ax3, var3; colormap = :balance, colorrange = (-1e-4, 1e-4))
heatmap!(ax4, var4; colormap = :balance, colorrange = (-1e-4, 1e-4))


title = @lift "Buoyancy, time = " * string(round(times[$n]./3600., digits=2)) * " hours"
Label(fig[1, 1:4], title, fontsize=24, tellwidth=false)

fig

frames = 1:length(times)

@info "Making an animation"

CairoMakie.record(fig, "geostrophic_adjustment_filtered_buoyancy_movie_offline.mp4", frames, framerate=24) do i
    n[] = i
end
# ![](geostrophic_adjustment_filtered_buoyancy_movie_offline.mp4)


# Then the tracer:
timeseries1 = FieldTimeSeries(filter_config.output_filename, "T")
timeseries2 = FieldTimeSeries(filter_config.output_filename, "T_Eulerian_filtered")
timeseries3 = FieldTimeSeries(filter_config.output_filename, "T_Lagrangian_filtered")
timeseries4 = FieldTimeSeries(filter_config.output_filename, "T_Lagrangian_filtered_at_mean")

times = timeseries1.times

set_theme!(Theme(fontsize = 20))
fig = Figure(size = (1300, 500))

axis_kwargs = (xlabel = "x",
               ylabel = "z",
               limits = ((-5000, 5000), (-100, 0)),
               aspect = AxisAspect(1))

ax1 = Axis(fig[2, 1]; title = "Raw", axis_kwargs...)
ax2 = Axis(fig[2, 2]; title = "Eulerian filtered", axis_kwargs...)
ax3 = Axis(fig[2, 3]; title = "Lagrangian filtered", axis_kwargs...)
ax4 = Axis(fig[2, 4]; title = "Lagrangian filtered \n at mean position", axis_kwargs...)


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


title = @lift "Tracer concentration, time = " * string(round(times[$n]./3600., digits=2)) * " hours"
Label(fig[1, 1:4], title, fontsize=24, tellwidth=false)

fig

frames = 1:length(times)

@info "Making an animation"

CairoMakie.record(fig, "geostrophic_adjustment_filtered_tracer_movie_offline.mp4", frames, framerate=24) do i
    n[] = i
end
# ![](geostrophic_adjustment_filtered_tracer_movie_offline.mp4)

# We see that the Eulerian filter smudges the tracer field as it is advected by the 
# inertial oscillations. The Lagrangian means directly calculated by this method are 
# identical to the raw fields for the tracer and buoyancy shown, as they are conservative  
# fields. However, when we remap to a mean reference position, we see the value of the
# Lagrangian filter in effectively removing the oscillations while preserving the tracer
# structures.

# Then the velocity into the page, v:
timeseries1 = FieldTimeSeries(filter_config.output_filename, "v")
timeseries2 = FieldTimeSeries(filter_config.output_filename, "v_Eulerian_filtered")
timeseries3 = FieldTimeSeries(filter_config.output_filename, "v_Lagrangian_filtered")
timeseries4 = FieldTimeSeries(filter_config.output_filename, "v_Lagrangian_filtered_at_mean")

times = timeseries1.times

set_theme!(Theme(fontsize = 20))
fig = Figure(size = (1300, 500))

axis_kwargs = (xlabel = "x",
               ylabel = "z",
               limits = ((-5000, 5000), (-100, 0)),
               aspect = AxisAspect(1))

ax1 = Axis(fig[2, 1]; title = "Raw", axis_kwargs...)
ax2 = Axis(fig[2, 2]; title = "Eulerian filtered", axis_kwargs...)
ax3 = Axis(fig[2, 3]; title = "Lagrangian filtered", axis_kwargs...)
ax4 = Axis(fig[2, 4]; title = "Lagrangian filtered \n at mean position", axis_kwargs...)


n = Observable(1)
Observable(1)

var1 = @lift timeseries1[$n]
var2 = @lift timeseries2[$n]
var3 = @lift timeseries3[$n]
var4 = @lift timeseries4[$n]

heatmap!(ax1, var1; colormap = :balance, colorrange = (-0.05, 0.05))
heatmap!(ax2, var2; colormap = :balance, colorrange = (-0.05, 0.05))
heatmap!(ax3, var3; colormap = :balance, colorrange = (-0.05, 0.05))
heatmap!(ax4, var4; colormap = :balance, colorrange = (-0.05, 0.05))


title = @lift "Velocity v, time = " * string(round(times[$n]./3600., digits=2)) * " hours"
Label(fig[1, 1:4], title, fontsize=24, tellwidth=false)

fig

frames = 1:length(times)

@info "Making an animation"

CairoMakie.record(fig, "geostrophic_adjustment_filtered_v_movie_offline.mp4", frames, framerate=24) do i
    n[] = i
end
# ![](geostrophic_adjustment_filtered_v_movie_offline.mp4)

# In this case, the Lagrangian filtered velocity fields differ from the raw fields, as expected, 
# since velocity is not a conservative field.

# We remove these files to keep things tidy, keep them for analysis if desired
rm(filename_stem * ".jld2")
rm(filter_config.output_filename)
rm(filter_config.output_filename[1:end-5] * ".nc")
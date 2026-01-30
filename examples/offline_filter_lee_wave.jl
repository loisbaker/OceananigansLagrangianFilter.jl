# # Steady lee waves with offline Lagrangian filtering

# We set up a two-dimensional hydrostatic simulation of a steady flow over a
# Gaussian bump, generating lee waves. We then apply an offline Lagrangian filter
# to the output of the simulation to extract the mean flow, removing the wave
# oscillations.

# The simulation is initialised with a constant and uniform stratification and horizontal background velocity

# In this example, the filtering is performed offline after the simulation. We let the simulation run to steady
# state, and run the filter on the steady data. This could also be performed online during the simulation, as 
# shown in the `lee_wave_online_filtering.jl` example.

# ## Run the simulation
# ### Install dependencies

using Oceananigans 
using Oceananigans.Units
using NCDatasets
using Printf
using Oceananigans.Grids: xnode, znode

# ### Model parameters

# Geometry
#Nx, Nz = 300, 400         # Higher res
Nx, Nz = 100, 100         # Lower res
H = 2kilometers           # Depth
L = 20kilometers          # Domain length

# Gaussian bump
h0 = 180meters            # Bump height
hill_width = 2kilometers  # Gaussian bump width 
x0 = -10kilometers        # Bump location
hill(x) = h0 * exp(-(x-x0)^2 / 2hill_width^2)
bottom(x) = - H + hill(x)

# Flow
f = 1e-4                  # Coriolis frequency [s⁻¹]
U = 0.2                   # Background flow speed [m/s]
Nsqr = 1e-6               # Buoyancy frequency squared [s⁻²]

κ = 0.1 # Diffusivity

filename_stem = "lee_wave";

# ### Define the grid

underlying_grid = RectilinearGrid(CPU(), size = (Nx, Nz), halo = (4, 4),
                                  x = (-L, L), z = (-H, 0),
                                  topology = (Periodic, Flat, Bounded))         

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom))

# ### Closures                     
κ = 0.1 
closure = ScalarDiffusivity(ν=κ, κ=κ)

# ### Forcing

# Steady geostrophic forcing through the v-momentum equation
A = f * U 
@inline steady_forcing(x, z, t, p) = p.A 
v_steady_forcing = Forcing(steady_forcing, parameters=(; A))

# Sponge layer to absorb waves at periodic boundaries
t_restore = 4hour
sponge_x1 = -4*L/5
sponge_x2 = 4*L/5
sponge_dx = L/5

sponge_mask(x, p) = 1.0 + 0.5*(tanh((x + p.sponge_x1)/p.sponge_dx) - tanh((x + p.sponge_x2)/p.sponge_dx))

# HydrostaticFreeSurfaceModel has issues with continuous forcings on GPU, so we use discrete_form=true
@inline function u_sponge_func(i, j, k, grid, clock, model_fields, p)
    timescale = p.t_restore
    x = xnode(i, j, k, grid, Face(), Center(), Center())
    u = @inbounds model_fields.u[i, j, k]

    return -1 / timescale * sponge_mask(x, p) * (u - p.U)
end

@inline function v_sponge_func(i, j, k, grid, clock, model_fields, p)
    timescale = p.t_restore
    x = xnode(i, j, k, grid, Center(), Face(), Center())
    v = @inbounds model_fields.v[i, j, k]

    return -1 / timescale * sponge_mask(x, p) * v 
end

@inline function b_sponge_func(i, j, k, grid, clock, model_fields, p)
    timescale = p.t_restore
    x = xnode(i, j, k, grid, Center(), Center(), Center())
    z = znode(i, j, k, grid, Center(), Center(), Center())
    b = @inbounds model_fields.b[i, j, k]

    return -1 / timescale * sponge_mask(x, p) * (b - p.Nsqr * z)
end

u_sponge_forcing = Forcing(u_sponge_func, discrete_form=true, parameters=(; U, t_restore, sponge_x1, sponge_x2, sponge_dx))
v_sponge_forcing = Forcing(v_sponge_func, discrete_form=true, parameters=(; t_restore, sponge_x1, sponge_x2, sponge_dx))
b_sponge_forcing = Forcing(b_sponge_func, discrete_form=true, parameters=(; Nsqr, t_restore, sponge_x1, sponge_x2, sponge_dx))


# ### Define the model
model = HydrostaticFreeSurfaceModel(grid; 
                                    coriolis = FPlane(f = f),
                                    buoyancy = BuoyancyTracer(),
                                    tracers = :b,
                                    tracer_advection = WENO(),
                                    momentum_advection = WENO(),
                                    forcing = (; v = (v_steady_forcing, v_sponge_forcing),
                                                u = u_sponge_forcing, b = b_sponge_forcing),                                                                                       
                                    closure=closure,
                                    )

# ### Initialise the buoyancy
bᵢ(x, z) = Nsqr * z
set!(model, u=U, b=bᵢ)

# ### Define the simulation
simulation = Simulation(model, Δt=20seconds, stop_time=15days)

# Set an adaptive timestep
conjure_time_step_wizard!(simulation, IterationInterval(20), cfl=0.2, max_Δt=2minutes)

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

add_callback!(simulation, print_progress, IterationInterval(200))


# ### Set up the output 
u, v, w = model.velocities
b = model.tracers.b

# Output a jld2 file for Lagrangian filtering
simulation.output_writers[:jld2fields] = JLD2Writer(
    model, (; u, w, b), filename = filename_stem * ".jld2", schedule=TimeInterval(1hour), overwrite_existing=true)

# ### Run simulation
@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)

# ## Perform Lagrangian filtering
# Now we set up and run the offline Lagrangian filter on the output of the above simulation.
# This could be performed in a different script (with appropriate import of Oceananigans.Units and CUDA if needed)

using OceananigansLagrangianFilter

# ### Set up the filter configuration
filter_config = OfflineFilterConfig(original_data_filename="lee_wave.jld2", # Where the original simulation output is
                                    output_filename = "lee_wave_offline_filtered.jld2", # Where to save the filtered output
                                    var_names_to_filter = ("b",), # Which variables to filter
                                    velocity_names = ("u","w"), # Velocities to use for remapping
                                    architecture = CPU(), # CPU() or GPU(), if GPU() make sure you have CUDA.jl installed and imported
                                    Δt = 1minutes, # Time step of filtering simulation
                                    T_start = 5days,
                                    T_out = 1hour, # How often to output filtered data
                                    N = 4, # Order of Butterworth filter
                                    freq_c = 1e-4/2, # Cut-off frequency of Butterworth filter
                                    compute_mean_velocities = true, # Whether to compute the mean velocities
                                    output_netcdf = true, # Whether to output filtered data to a netcdf file in addition to .jld2
                                    delete_intermediate_files = true, # Delete the individual output of the forward and backward passes
                                    compute_Eulerian_filter = true) # Whether to compute the Eulerian filter for comparison


# ### Run the offline Lagrangian filter
run_offline_Lagrangian_filter(filter_config)

# ### Visualisation

using CairoMakie 

# Now we animate the filterieng results for the horizontal velocity u and buoyancy b:
timeseries1 = FieldTimeSeries(filter_config.output_filename, "u")
timeseries2 = FieldTimeSeries(filter_config.output_filename, "u_Eulerian_filtered")
timeseries3 = FieldTimeSeries(filter_config.output_filename, "u_Lagrangian_filtered")
timeseries4 = FieldTimeSeries(filter_config.output_filename, "u_Lagrangian_filtered_at_mean")

b_timeseries1 = FieldTimeSeries(filter_config.output_filename, "b")
b_timeseries2 = FieldTimeSeries(filter_config.output_filename, "b_Eulerian_filtered")
b_timeseries3 = FieldTimeSeries(filter_config.output_filename, "b_Lagrangian_filtered")
b_timeseries4 = FieldTimeSeries(filter_config.output_filename, "b_Lagrangian_filtered_at_mean")

times = timeseries1.times
bottom_height = vec(timeseries1.grid.immersed_boundary.bottom_height)
Nx = timeseries1.grid.underlying_grid.Nx
x = Array(timeseries1.grid.underlying_grid.xᶜᵃᵃ[1:Nx])

set_theme!(Theme(fontsize = 25))
fig = Figure(size = (2000, 950))

axis_kwargs = (xlabel = "x",
               ylabel = "z",
               limits = ((-20000, 10000), (-2000, 0)),
               aspect = AxisAspect(2.5))

ax1 = Axis(fig[2, 1]; title = "Raw", axis_kwargs...)
ax2 = Axis(fig[2, 2]; title = "Eulerian filtered", axis_kwargs...)
ax3 = Axis(fig[3, 1]; title = "Lagrangian filtered", axis_kwargs...)
ax4 = Axis(fig[3, 2]; title = "Lagrangian filtered at mean", axis_kwargs...)

n = Observable(1)
Observable(1)

var1 = @lift timeseries1[$n]
var2 = @lift timeseries2[$n]
var3 = @lift timeseries3[$n]
var4 = @lift timeseries4[$n]
b_var1 = @lift b_timeseries1[$n]
b_var2 = @lift b_timeseries2[$n]
b_var3 = @lift b_timeseries3[$n]
b_var4 = @lift b_timeseries4[$n]


hm_1 = heatmap!(ax1, var1; colormap = :balance, colorrange = (0.1, 0.3))
hm_2 = heatmap!(ax2, var2; colormap = :balance, colorrange = (0.1, 0.3))
hm_3 = heatmap!(ax3, var3; colormap = :balance, colorrange = (0.19, 0.21))
hm_4 = heatmap!(ax4, var4; colormap = :balance, colorrange = (0.19, 0.21))

contour!(ax1, b_var1; levels = 40, color = :black, linewidth = 0.5)
contour!(ax2, b_var2; levels = 40, color = :black, linewidth = 0.5)
contour!(ax3, b_var3; levels = 40, color = :black, linewidth = 0.5)
contour!(ax4, b_var4; levels = 40, color = :black, linewidth = 0.5)

Colorbar(fig[2, 3], hm_1, label = "m s⁻¹")
Colorbar(fig[3, 3], hm_3, label = "m s⁻¹")

# Plot topography
floor_level = fill(-2000, length(x)) 
band!(ax1, x, floor_level, bottom_height, color=:black)
band!(ax2, x, floor_level, bottom_height, color=:black)
band!(ax3, x, floor_level, bottom_height, color=:black)
band!(ax4, x, floor_level, bottom_height, color=:black)


title = @lift "Horizontal velocity, time = " * string(round(times[$n]./3600/24, digits=2)) * " days"
Label(fig[1, 1:2], title, fontsize=30, tellwidth=false)


fig

frames = 1:length(times)

@info "Making an animation"

CairoMakie.record(fig, "lee_wave_filtered_u_movie_offline.mp4", frames, framerate=24) do i
    n[] = i
end
# ![](lee_wave_filtered_u_movie_offline.mp4)

# The Eulerian filtered output looks very similar to the raw output, since the lee waves are steady. 
# They're steady because they've been Doppler-shifted by the mean flow. When we use the Lagrangian 
# filter instead, we see that the lee waves are removed as they are high frequency in the frame of the 
# flow. Note that the colour range shown is an order of magnitude smaller for the Lagrangian filtered 
# fields, allowing us to see the wave impact on the mean flow. The contours of buoyancy demonstrate 
# the difference between the displaced Lagrangian filtered field (bottom left) and the Lagrangian 
#filtered field remapped to the mean position (bottom right). This removes the wave-like displacements
# from the isopycnals. The interpolation stage to calculate this field is imperfect near the immersed 
# boundary (it's not guaranteed that every spatial location has a corresponding Lagrangian mean defined
# at the mean position), so it has been masked.


# We remove these files to keep things tidy, keep them for analysis if desired
rm(filename_stem * ".jld2")
rm(filter_config.output_filename)
rm(filter_config.output_filename[1:end-5] * ".nc")
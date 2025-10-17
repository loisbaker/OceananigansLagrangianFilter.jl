using OceananigansLagrangianFilter
using OceananigansLagrangianFilter.Utils # contains some utility functions
using Oceananigans.Units
using Oceananigans.TurbulenceClosures
using CairoMakie 
using NCDatasets
using Printf

# We set up a geostrophic adjustment problem similar to Blumen (2000), JPO
# in a domain that is horizontally periodic, and Lagrangian filter it as we go. 
# Credit to Tom Cummings for work on this example. 

# Model parameters
Nx = 300
Nz = 80
f = 1e-4                # Coriolis frequency [s⁻¹]
L_front = 10kilometers  # Initial front width [m]
aspect_ratio = 100      # L_front/H
Ro = 0.1                # Rossby number (defines M^2)

# Derived parameters
H = L_front/aspect_ratio  # Depth
M² = (Ro^2*f^2*L_front)/H # Horizontal buoyancy gradient
Δb = M²*L_front # Buoyancy difference across the front
κh = 1e-4 # Horizontal diffusivity
κv = 1e-4 # Vertical diffusivity

filename_stem = "geostrophic_adjustment"

grid = RectilinearGrid(CPU(),size = (Nx, Nz), 
                       x = (-L_front/2, L_front/2),
                       z = (-H, 0),
                       topology = (Periodic, Flat, Bounded))

# Define tracers
tracers = (:b,:T)

# Define any forcing
forcing = NamedTuple()

# Define filter configuration
filter_config = OnlineFilterConfig( grid = grid,
                                    output_filename = filename_stem * ".jld2",
                                    var_names_to_filter = ("b","T"), 
                                    velocity_names = ("u","w"),
                                    N = 2,
                                    freq_c = f/2)

# Create the filtered variables - these will be tracers in the model
filtered_vars = create_filtered_vars(filter_config)

# Add to the existing tracers
tracers = (filtered_vars..., tracers...)

# Create forcing for these filtered variables
filter_forcing = create_forcing(filtered_vars, filter_config)

# Add to the existing forcing
forcing = merge(forcing, filter_forcing);

# No real need for a closure here, but we include one for demo purposes

# Helper to set filtered variable closures to zero
zero_filtered_var_closure = zero_closure_for_filtered_vars(filter_config)

# Unless we set filtered variable closures to zero, the default closure will apply to all tracers
horizontal_closure = HorizontalScalarDiffusivity(ν=1e-6, κ=merge((T=1e-6, b= 1e-6),zero_filtered_var_closure) )
vertical_closure = VerticalScalarDiffusivity(ν=1e-6 , κ=merge((T=1e-6, b= 1e-6),zero_filtered_var_closure) )
closure = (horizontal_closure, vertical_closure)

# Define the model
model =  NonhydrostaticModel(; grid,
                coriolis = FPlane(f = f),
                buoyancy = BuoyancyTracer(),
                tracers = tracers,
                forcing = forcing,
                advection = WENO(),
                closure = closure)

# Initialise the buoyancy and tracer (velocities start at rest by default)
bᵢ(x, z) = Δb*sin(2*pi/L_front * x)
Tᵢ(x, z) = exp(-(x/(L_front/50)).^2)
set!(model, b= bᵢ, T= Tᵢ) 

# Set appropriate initial conditions for the filtered variables based on the actual variables
initialise_filtered_vars_from_model(model, filter_config)

# Define the simulation
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

# Set up the output 

# Create filtered outputs
outputs = create_output_fields(model, filter_config)

# Add in original variables if needed
outputs["b"] = model.tracers.b
outputs["T"] = model.tracers.T
outputs["u"] = model.velocities.u
outputs["v"] = model.velocities.v
outputs["w"] = model.velocities.w

# Outout a jld2 file
simulation.output_writers[:jld2fields] = JLD2Writer(
    model, outputs, filename=filter_config.output_filename, schedule=TimeInterval(1hour), overwrite_existing=true)


# Could also output a netcdf file
rm(filename_stem * ".nc",force=true)
simulation.output_writers[:ncfields] = NetCDFWriter(
    model, outputs, filename=filename_stem * ".nc", schedule=TimeInterval(1hour), overwrite_existing=true)
   
# Run the simulation
@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)

# Option to regrid to mean position
if filter_config.map_to_mean
    regrid_to_mean_position!(filter_config)
end

# Option to calculate Eulerian filter too
compute_Eulerian_filter!(filter_config)


# Option to add a shifted time coordinate
compute_time_shift!(filter_config)


# Option to output final netcdf
jld2_to_netcdf(filename_stem * ".jld2", filename_stem * ".nc")


# Animate the results, buoyancy first:
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


title = @lift "Buoyancy, time = " * string(round(times[$n], digits=2))
Label(fig[1, 1:4], title, fontsize=24, tellwidth=false)

fig

frames = 1:length(times)

@info "Making an animation"

CairoMakie.record(fig, "geostrophic_adjustment_filtered_buoyancy_movie_online.mp4", frames, framerate=24) do i
    n[] = i
end


# Then the tracer:
# Animate
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


title = @lift "Tracer concentration, time = " * string(round(times[$n], digits=2))
Label(fig[1, 1:4], title, fontsize=24, tellwidth=false)

fig

frames = 1:length(times)

@info "Making an animation"

CairoMakie.record(fig, "geostrophic_adjustment_filtered_tracer_movie_online.mp4", frames, framerate=24) do i
    n[] = i
end
# # Geostrophic adjustment with online Lagrangian filtering

# We set up a geostrophic adjustment problem similar to Blumen (2000), *JPO*
# in a domain that is horizontally periodic. 

# An initially unbalanced two-dimensional front oscillates with the inertial 
# frequency around a state of geostrophic balance, and we illustrate that we 
# can remove the oscillations to find the mean state. Thanks to Tom Cummings 
# for work on this example. 

# In this example, the filtering is performed online during the simulation.

# ### Load dependencies
using OceananigansLagrangianFilter # Gives access to all Oceananigans functions too
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

filename_stem = "geostrophic_adjustment_online_filtered";

# ### Grid
grid = RectilinearGrid(CPU(),size = (Nx, Nz), 
                       x = (-L_front/2, L_front/2),
                       z = (-H, 0),
                       topology = (Periodic, Flat, Bounded))

# ### Define model tracers
tracers = (:b,:T);

# ### Define model forcing
forcing = NamedTuple();

# ### Define filter configuration
filter_config = OnlineFilterConfig( grid = grid,
                                    output_filename = filename_stem * ".jld2",
                                    var_names_to_filter = ("b","T"), 
                                    velocity_names = ("u","w","v"),
                                    N = 2,
                                    freq_c = f/2)

# ### Create the filtered variables - these will be tracers in the model
filtered_vars = create_filtered_vars(filter_config)

# ### Add to the existing tracers
tracers = (filtered_vars..., tracers...)

# ### Create forcing for these filtered variables
filter_forcing = create_forcing(filtered_vars, filter_config);

# ### Add these to the existing forcing
forcing = merge(forcing, filter_forcing);

# ### Define closures
# If the model uses a closure, we use a helper function to set filtered variable closures to zero (unless we set filtered variable closures to zero, the default closure will apply to all tracers).
zero_filtered_var_closure = zero_closure_for_filtered_vars(filter_config)
horizontal_closure = HorizontalScalarDiffusivity(ν=κh, κ=merge((T=κh, b= κh),zero_filtered_var_closure) )
vertical_closure = VerticalScalarDiffusivity(ν=κv , κ=merge((T=κv, b= κv),zero_filtered_var_closure) )
closure = (horizontal_closure, vertical_closure);
nothing #hide

# ### Define the model
model =  NonhydrostaticModel(; grid,
                coriolis = FPlane(f = f),
                buoyancy = BuoyancyTracer(),
                tracers = tracers,
                forcing = forcing,
                advection = WENO(),
                closure = closure)

# ### Initialise tracers
# Model buoyancy and tracers
bᵢ(x, z) = Δb*sin(2*pi/L_front * x)
Tᵢ(x, z) = exp(-(x/(L_front/50)).^2)
set!(model, b= bᵢ, T= Tᵢ) 

# Set appropriate initial conditions for the filtered variables based on the actual variables
initialise_filtered_vars_from_model(model, filter_config)

# ### Define the simulation
simulation = Simulation(model, Δt=20minutes, stop_time=3days)

# ### Set an adaptive timestep
conjure_time_step_wizard!(simulation, IterationInterval(20), cfl=0.2, max_Δt=20minutes)

# ### Add a progress callback

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

# ### Set up the outputs
# Create filtered outputs:
outputs = create_output_fields(model, filter_config);

# Add in original variables if needed:
outputs["b"] = model.tracers.b;
outputs["T"] = model.tracers.T;
outputs["u"] = model.velocities.u;
outputs["v"] = model.velocities.v;
outputs["w"] = model.velocities.w;

# Output a .jld2 file:
simulation.output_writers[:jld2fields] = JLD2Writer(
    model, outputs, filename=filter_config.output_filename, schedule=TimeInterval(1hour), overwrite_existing=true)

# ### Run the simulation
@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)

# ### Option to regrid to mean position
if filter_config.map_to_mean
    regrid_to_mean_position!(filter_config)
end
nothing #hide

# ### Option to calculate Eulerian filter too
compute_Eulerian_filter!(filter_config);
nothing #hide

# ### Option to add a shifted time coordinate
compute_time_shift!(filter_config)


# ### Option to output a final NetCDF file
jld2_to_netcdf(filename_stem * ".jld2", filename_stem * ".nc")


# ### Animate the results, buoyancy first:

using CairoMakie 

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

CairoMakie.record(fig, "geostrophic_adjustment_filtered_buoyancy_movie_online.mp4", frames, framerate=24) do i
    n[] = i
end
# ![](geostrophic_adjustment_filtered_buoyancy_movie_online.mp4)


# ### Then plot the tracer concentration:

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


title = @lift "Tracer concentration, time = " * string(round(times[$n]./3600., digits=2)) * " hours"
Label(fig[1, 1:4], title, fontsize=24, tellwidth=false)

fig

frames = 1:length(times)

@info "Making an animation"

CairoMakie.record(fig, "geostrophic_adjustment_filtered_tracer_movie_online.mp4", frames, framerate=24) do i
    n[] = i
end
nothing #hide

# ![](geostrophic_adjustment_filtered_tracer_movie_online.mp4)


# We see that the Eulerian filter smudges the tracer field as it is advected by the 
# inertial oscillations. The Lagrangian means directly calculated by this method are 
# identical to the raw fields for the tracer and buoyancy shown, as they are conservative  
# fields. However, when we remap to a mean reference position, we see the value of the
# Lagrangian filter in effectively removing the oscillations while preserving the tracer
# structures. In comparison to the offline filtering example, the online filter does a 
# slightly worse job removing the oscillations in the 'Lagrangian filtered at mean' field,
# since the filter is a more optimal low-pass.

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

CairoMakie.record(fig, "geostrophic_adjustment_filtered_v_movie_online.mp4", frames, framerate=24) do i
    n[] = i
end
# ![](geostrophic_adjustment_filtered_v_movie_online.mp4)

# In this case, the Lagrangian filtered velocity fields differ from the raw fields, as expected, 
# since velocity is not a conservative field.

# We remove these files to keep things tidy, keep them for analysis if desired
rm(filename_stem * ".jld2")
rm(filename_stem * ".nc")
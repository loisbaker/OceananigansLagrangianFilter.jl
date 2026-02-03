```@meta
EditURL = "../../../examples/online_filter_geostrophic_adjustment.jl"
```

# Geostrophic adjustment with online Lagrangian filtering

We set up a geostrophic adjustment problem similar to Blumen (2000), *JPO*
in a domain that is horizontally periodic.

An initially unbalanced two-dimensional front oscillates with the inertial
frequency around a state of geostrophic balance, and we illustrate that we
can remove the oscillations to find the mean state. Thanks to Tom Cummings
for work on this example.

In this example, the filtering is performed online during the simulation.

### Load dependencies

````julia
using OceananigansLagrangianFilter # Gives access to all Oceananigans functions too
using Oceananigans.Units
using NCDatasets
using Printf
````

### Model parameters

````julia
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
````

### Grid

````julia
grid = RectilinearGrid(CPU(),size = (Nx, Nz),
                       x = (-L_front/2, L_front/2),
                       z = (-H, 0),
                       topology = (Periodic, Flat, Bounded))
````

````
400×1×80 RectilinearGrid{Float64, Periodic, Flat, Bounded} on CPU with 3×0×3 halo
├── Periodic x ∈ [-5000.0, 5000.0) regularly spaced with Δx=25.0
├── Flat y                         
└── Bounded  z ∈ [-100.0, 0.0]     regularly spaced with Δz=1.25
````

### Define model tracers

````julia
tracers = (:b,:T);
````

### Define model forcing

````julia
forcing = NamedTuple();
````

### Define filter configuration

````julia
filter_config = OnlineFilterConfig( grid = grid,
                                    output_filename = filename_stem * ".jld2",
                                    var_names_to_filter = ("b","T"),
                                    velocity_names = ("u","w","v"),
                                    N = 2,
                                    freq_c = f/4)
````

````
OnlineFilterConfig(400×1×80 RectilinearGrid{Float64, Periodic, Flat, Bounded} on CPU with 3×0×3 halo
├── Periodic x ∈ [-5000.0, 5000.0) regularly spaced with Δx=25.0
├── Flat y                         
└── Bounded  z ∈ [-100.0, 0.0]     regularly spaced with Δz=1.25, "geostrophic_adjustment_online_filtered.jld2", ("b", "T"), ("u", "w", "v"), (a1 = 7.10533784274036e-21, b1 = -3.535533905932738e-5, c1 = 1.767766952966369e-5, d1 = -1.767766952966369e-5, N_coeffs = 1), true, true, 5, "online", "")
````

### Create the filtered variables - these will be tracers in the model

````julia
filtered_vars = create_filtered_vars(filter_config)
````

````
(:b_C1, :T_C1, :xi_u_C1, :xi_w_C1, :xi_v_C1, :b_S1, :T_S1, :xi_u_S1, :xi_w_S1, :xi_v_S1)
````

### Add to the existing tracers

````julia
tracers = (filtered_vars..., tracers...)
````

````
(:b_C1, :T_C1, :xi_u_C1, :xi_w_C1, :xi_v_C1, :b_S1, :T_S1, :xi_u_S1, :xi_w_S1, :xi_v_S1, :b, :T)
````

### Create forcing for these filtered variables

````julia
filter_forcing = create_forcing(filtered_vars, filter_config);
````

### Add these to the existing forcing

````julia
forcing = merge(forcing, filter_forcing);
````

### Define closures
If the model uses a closure, we use a helper function to set filtered variable closures to zero (unless we set filtered variable closures to zero, the default closure will apply to all tracers).

````julia
zero_filtered_var_closure = zero_closure_for_filtered_vars(filter_config)
horizontal_closure = HorizontalScalarDiffusivity(ν=κh, κ=merge((T=κh, b= κh),zero_filtered_var_closure) )
vertical_closure = VerticalScalarDiffusivity(ν=κv , κ=merge((T=κv, b= κv),zero_filtered_var_closure) )
closure = (horizontal_closure, vertical_closure);
````

### Define the model

````julia
model =  NonhydrostaticModel(grid;
                coriolis = FPlane(f = f),
                buoyancy = BuoyancyTracer(),
                tracers = tracers,
                forcing = forcing,
                advection = WENO(),
                closure = closure)
````

````
NonhydrostaticModel{CPU, RectilinearGrid}(time = 0 seconds, iteration = 0)
├── grid: 400×1×80 RectilinearGrid{Float64, Periodic, Flat, Bounded} on CPU with 3×0×3 halo
├── timestepper: RungeKutta3TimeStepper
├── advection scheme: WENO{3, Float64, Float32}(order=5)
├── tracers: (b_C1, T_C1, xi_u_C1, xi_w_C1, xi_v_C1, b_S1, T_S1, xi_u_S1, xi_w_S1, xi_v_S1, b, T)
├── closure: Tuple with 2 closures:
│   ├── HorizontalScalarDiffusivity{ExplicitTimeDiscretization}(ν=1.0e-6, κ=(b_C1=0.0, T_C1=0.0, xi_u_C1=0.0, xi_w_C1=0.0, xi_v_C1=0.0, b_S1=0.0, T_S1=0.0, xi_u_S1=0.0, xi_w_S1=0.0, xi_v_S1=0.0, b=1.0e-6, T=1.0e-6))
│   └── VerticalScalarDiffusivity{ExplicitTimeDiscretization}(ν=1.0e-6, κ=(b_C1=0.0, T_C1=0.0, xi_u_C1=0.0, xi_w_C1=0.0, xi_v_C1=0.0, b_S1=0.0, T_S1=0.0, xi_u_S1=0.0, xi_w_S1=0.0, xi_v_S1=0.0, b=1.0e-6, T=1.0e-6))
├── buoyancy: BuoyancyTracer with ĝ = NegativeZDirection()
└── coriolis: FPlane{Float64}(f=0.0001)
````

### Initialise tracers
Model buoyancy and tracers

````julia
bᵢ(x, z) = Δb*sin(2*pi/L_front * x)
Tᵢ(x, z) = exp(-(x/(L_front/50)).^2)
set!(model, b= bᵢ, T= Tᵢ)
````

Set appropriate initial conditions for the filtered variables based on the actual variables

````julia
initialise_filtered_vars_from_model(model, filter_config)
````

### Define the simulation

````julia
simulation = Simulation(model, Δt=20minutes, stop_time=3days)
````

````
Simulation of NonhydrostaticModel{CPU, RectilinearGrid}(time = 0 seconds, iteration = 0)
├── Next time step: 20 minutes
├── run_wall_time: 0 seconds
├── run_wall_time / iteration: NaN days
├── stop_time: 3 days
├── stop_iteration: Inf
├── wall_time_limit: Inf
├── minimum_relative_step: 0.0
├── callbacks: OrderedDict with 4 entries:
│   ├── stop_time_exceeded => Callback of stop_time_exceeded on IterationInterval(1)
│   ├── stop_iteration_exceeded => Callback of stop_iteration_exceeded on IterationInterval(1)
│   ├── wall_time_limit_exceeded => Callback of wall_time_limit_exceeded on IterationInterval(1)
│   └── nan_checker => Callback of NaNChecker for u on IterationInterval(100)
└── output_writers: OrderedDict with no entries
````

### Set an adaptive timestep

````julia
conjure_time_step_wizard!(simulation, IterationInterval(20), cfl=0.2, max_Δt=20minutes)
````

### Add a progress callback

````julia
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
````

### Set up the outputs
Create filtered outputs:

````julia
outputs = create_output_fields(model, filter_config);
````

Add in original variables if needed:

````julia
outputs["b"] = model.tracers.b;
outputs["T"] = model.tracers.T;
outputs["u"] = model.velocities.u;
outputs["v"] = model.velocities.v;
outputs["w"] = model.velocities.w;
````

Output a .jld2 file:

````julia
simulation.output_writers[:jld2fields] = JLD2Writer(
    model, outputs, filename=filter_config.output_filename, schedule=TimeInterval(1hour), overwrite_existing=true)
````

````
JLD2Writer scheduled on TimeInterval(1 hour):
├── filepath: geostrophic_adjustment_online_filtered.jld2
├── 13 outputs: (xi_u, T, b, b_Lagrangian_filtered, xi_v, xi_w, v, w, v_Lagrangian_filtered, w_Lagrangian_filtered, T_Lagrangian_filtered, u, u_Lagrangian_filtered)
├── array_type: Array{Float32}
├── including: [:grid, :coriolis, :buoyancy, :closure]
├── file_splitting: NoFileSplitting
└── file size: 0 bytes (file not yet created)
````

### Run the simulation

````julia
@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)
````

````
[ Info: Running the simulation...
[ Info: Initializing simulation...
[00.00%] i: 0, t: 0 seconds, wall time: 19.234 seconds, max(u): (0.000e+00, 0.000e+00, 0.000e+00) m/s, next Δt: 20 minutes
[ Info:     ... simulation initialization complete (4.321 minutes)
[ Info: Executing initial time step...
[ Info:     ... initial time step complete (16.058 seconds).
[15.28%] i: 50, t: 11 hours, wall time: 6.640 minutes, max(u): (2.205e-02, 5.272e-02, 3.914e-04) m/s, next Δt: 6.183 minutes
[19.31%] i: 100, t: 13.906 hours, wall time: 2.080 minutes, max(u): (2.981e-02, 2.366e-02, 5.021e-04) m/s, next Δt: 2.795 minutes
[22.61%] i: 150, t: 16.282 hours, wall time: 2.156 minutes, max(u): (1.448e-02, 3.717e-03, 1.989e-04) m/s, next Δt: 3.382 minutes
[26.86%] i: 200, t: 19.341 hours, wall time: 2.094 minutes, max(u): (1.885e-02, 6.552e-03, 3.549e-04) m/s, next Δt: 4.420 minutes
[30.87%] i: 250, t: 22.226 hours, wall time: 2.165 minutes, max(u): (3.049e-02, 3.460e-02, 5.432e-04) m/s, next Δt: 2.712 minutes
[34.20%] i: 300, t: 1.026 days, wall time: 2.078 minutes, max(u): (1.763e-02, 5.629e-02, 3.552e-04) m/s, next Δt: 3.422 minutes
[38.36%] i: 350, t: 1.151 days, wall time: 2.135 minutes, max(u): (1.565e-02, 5.854e-02, 3.038e-04) m/s, next Δt: 4.141 minutes
[42.41%] i: 400, t: 1.272 days, wall time: 2.061 minutes, max(u): (3.055e-02, 3.396e-02, 6.365e-04) m/s, next Δt: 2.727 minutes
[45.61%] i: 450, t: 1.368 days, wall time: 2.099 minutes, max(u): (2.307e-02, 1.076e-02, 3.970e-04) m/s, next Δt: 3.076 minutes
[49.47%] i: 500, t: 1.484 days, wall time: 2.055 minutes, max(u): (9.182e-03, 3.661e-03, 2.610e-04) m/s, next Δt: 4.094 minutes
[53.85%] i: 550, t: 1.616 days, wall time: 2.116 minutes, max(u): (2.978e-02, 2.383e-02, 6.374e-04) m/s, next Δt: 2.965 minutes
[57.01%] i: 600, t: 1.710 days, wall time: 2.087 minutes, max(u): (2.578e-02, 4.763e-02, 5.998e-04) m/s, next Δt: 3.092 minutes
[60.82%] i: 650, t: 1.825 days, wall time: 2.125 minutes, max(u): (5.830e-03, 6.150e-02, 1.671e-04) m/s, next Δt: 3.741 minutes
[65.28%] i: 700, t: 1.958 days, wall time: 2.077 minutes, max(u): (2.783e-02, 4.469e-02, 7.864e-04) m/s, next Δt: 2.994 minutes
[68.51%] i: 750, t: 2.055 days, wall time: 2.102 minutes, max(u): (2.840e-02, 2.019e-02, 7.026e-04) m/s, next Δt: 2.800 minutes
[72.01%] i: 800, t: 2.160 days, wall time: 2.074 minutes, max(u): (1.031e-02, 6.012e-03, 1.518e-04) m/s, next Δt: 3.727 minutes
[76.49%] i: 850, t: 2.295 days, wall time: 2.134 minutes, max(u): (2.391e-02, 1.160e-02, 7.266e-04) m/s, next Δt: 4.419 minutes
[80.12%] i: 900, t: 2.404 days, wall time: 2.130 minutes, max(u): (2.944e-02, 3.850e-02, 8.380e-04) m/s, next Δt: 2.831 minutes
[83.57%] i: 950, t: 2.507 days, wall time: 2.208 minutes, max(u): (1.288e-02, 5.856e-02, 5.939e-04) m/s, next Δt: 3.425 minutes
[87.88%] i: 1000, t: 2.637 days, wall time: 2.168 minutes, max(u): (2.045e-02, 5.464e-02, 7.651e-04) m/s, next Δt: 4.074 minutes
[91.73%] i: 1050, t: 2.752 days, wall time: 2.186 minutes, max(u): (3.010e-02, 2.809e-02, 1.100e-03) m/s, next Δt: 2.759 minutes
[95.10%] i: 1100, t: 2.853 days, wall time: 2.123 minutes, max(u): (1.714e-02, 8.333e-03, 4.073e-04) m/s, next Δt: 3.453 minutes
[99.29%] i: 1150, t: 2.979 days, wall time: 2.171 minutes, max(u): (1.560e-02, 8.750e-03, 8.093e-04) m/s, next Δt: 4.178 minutes
[ Info: Simulation is stopping after running for 53.638 minutes.
[ Info: Simulation time 3 days equals or exceeds stop time 3 days.
[ Info: Simulation completed in 53.690 minutes

````

### Option to regrid to mean position

````julia
if filter_config.map_to_mean
    regrid_to_mean_position!(filter_config)
end
````

````
[ Info: Assuming velocities normal to z boundaries are zero (open boundaries not yet supported for regridding)
[ Info: Wrote regridded data to new variables with _at_mean suffix in file geostrophic_adjustment_online_filtered.jld2

````

### Option to calculate Eulerian filter too

````julia
compute_Eulerian_filter!(filter_config);
````

````
[ Info: Computing Eulerian filter for variable b
[ Info: Computing Eulerian filter for variable T
[ Info: Computing Eulerian filter for variable u
[ Info: Computing Eulerian filter for variable w
[ Info: Computing Eulerian filter for variable v

````

### Option to add a shifted time coordinate

````julia
compute_time_shift!(filter_config)
````

````
[ Info: Wrote time shift data to new group timeseries/t_shifted in file geostrophic_adjustment_online_filtered.jld2

````

### Option to output a final NetCDF file

````julia
jld2_to_netcdf(filename_stem * ".jld2", filename_stem * ".nc")
````

````
[ Info: Wrote NetCDF file to geostrophic_adjustment_online_filtered.nc

````

### Animate the results, buoyancy first:

````julia
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
````

````
"geostrophic_adjustment_filtered_buoyancy_movie_online.mp4"
````

![](geostrophic_adjustment_filtered_buoyancy_movie_online.mp4)

### Then plot the tracer concentration:

````julia
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
````

````
[ Info: Making an animation

````

![](geostrophic_adjustment_filtered_tracer_movie_online.mp4)

We see that the Eulerian filter smudges the tracer field as it is advected by the
inertial oscillations. The Lagrangian means directly calculated by this method are
identical to the raw fields for the tracer and buoyancy shown, as they are conservative
fields. However, when we remap to a mean reference position, we see the value of the
Lagrangian filter in effectively removing the oscillations while preserving the tracer
structures. In comparison to the offline filtering example, the online filter does a
slightly worse job removing the oscillations in the 'Lagrangian filtered at mean' field,
since the filter is a more optimal low-pass.

Then the velocity into the page, v:

````julia
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
````

````
"geostrophic_adjustment_filtered_v_movie_online.mp4"
````

![](geostrophic_adjustment_filtered_v_movie_online.mp4)

In this case, the Lagrangian filtered velocity fields differ from the raw fields, as expected,
since velocity is not a conservative field.

````julia
# We remove these files to keep things tidy, keep them for analysis if desired
rm(filename_stem * ".jld2")
rm(filename_stem * ".nc")
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


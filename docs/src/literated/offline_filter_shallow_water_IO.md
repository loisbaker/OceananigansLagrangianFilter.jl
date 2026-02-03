```@meta
EditURL = "../../../examples/offline_filter_shallow_water_IO.jl"
```

# Shallow water intertial oscillations with offline Lagrangian filtering

This example demonstrates how to perform offline filtering on a shallow water simulation to
remove the effect of inertial oscillations on a tracer field.

We could also filter a shallow water simulation online, but would have to use the
``VectorInvariantFormulation`` in order to have direct access to the model velocities.
This example uses the ``ConservativeFormulation`` instead, and filtering is performed
offline using the saved velocities after they have been calculated from ``uh`` and ``vh``.

## Run the simulation
### Install dependencies

````julia
using Oceananigans
using Printf
using NCDatasets

filename_stem = "SW_IO_with_tracer";
````

### Define the grid

````julia
grid = RectilinearGrid(CPU(), size = (50, 50),
                       x = (0, 2*pi),
                       y = (0, 2*pi),
                       topology = (Periodic, Periodic, Flat))
````

````
50×50×1 RectilinearGrid{Float64, Periodic, Periodic, Flat} on CPU with 3×3×0 halo
├── Periodic x ∈ [4.03717e-17, 6.28319) regularly spaced with Δx=0.125664
├── Periodic y ∈ [4.03717e-17, 6.28319) regularly spaced with Δy=0.125664
└── Flat z                              
````

### Set parameters
Building a `ShallowWaterModel`. We non-dimensionalise as in Kafiabad & Vanneste 2023.

````julia
Fr = 0.1 # Froude number
Ro = 1 # fRossby number

gravitational_acceleration = 1/Fr^2
coriolis = FPlane(f=1/Ro)
````

````
FPlane{Float64}(f=1.0)
````

### Define the model

````julia
model = ShallowWaterModel(grid; coriolis, gravitational_acceleration,
                            timestepper = :RungeKutta3,
                            tracers= (:T,),
                            momentum_advection = WENO())
````

````
ShallowWaterModel{CPU, Float64}(time = 0 seconds, iteration = 0) 
├── grid: 50×50×1 RectilinearGrid{Float64, Periodic, Periodic, Flat} on CPU with 3×3×0 halo
├── timestepper: RungeKutta3TimeStepper
├── advection scheme: 
│   ├── momentum: WENO{3, Float64, Float32}(order=5)
│   ├── mass: WENO{3, Float64, Float32}(order=5)
│   └── T: WENO{3, Float64, Float32}(order=5)
├── tracers: (:T,)
└── coriolis: FPlane{Float64}
````

### Initial conditions

Velocity and height initial conditions - uniform velocity perturbation, initial height is 1 (unperturbed)

````julia
displacement = 2*pi/10
u_i = displacement/Ro
h_i = 1
uh_i = u_i*h_i;
````

Initialise a tracer as a blob in the middle of the domain

````julia
width = 2*pi/15
T_i(x, y) = exp(-(((x - pi)^2 + (y - pi)^2)/width).^2)

set!(model, uh = uh_i, h= h_i, T = T_i )

uh, vh, h = model.solution

u = Field(uh / h)
v = Field(vh / h)
T = model.tracers.T
````

````
50×50×1 Field{Center, Center, Center} on RectilinearGrid on CPU
├── grid: 50×50×1 RectilinearGrid{Float64, Periodic, Periodic, Flat} on CPU with 3×3×0 halo
├── boundary conditions: FieldBoundaryConditions
│   └── west: Periodic, east: Periodic, south: Periodic, north: Periodic, bottom: Nothing, top: Nothing, immersed: Nothing
└── data: 56×56×1 OffsetArray(::Array{Float64, 3}, -2:53, -2:53, 1:1) with eltype Float64 with indices -2:53×-2:53×1:1
    └── max=0.999645, min=0.0, mean=0.0295409
````

### Simulation

````julia
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
````

````
Callback of progress on IterationInterval(100)
````

### Set up the outputs
Save velocities and tracer for Lagrangian filtering

````julia
simulation.output_writers[:fields_jld2] = JLD2Writer(model, (; u,v,T),
                                                        filename = filename_stem * ".jld2",
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)
````

````
JLD2Writer scheduled on TimeInterval(100 ms):
├── filepath: SW_IO_with_tracer.jld2
├── 3 outputs: (u, v, T)
├── array_type: Array{Float32}
├── including: [:grid, :coriolis, :closure]
├── file_splitting: NoFileSplitting
└── file size: 0 bytes (file not yet created)
````

### And finally run the simulation.

````julia
run!(simulation)
````

````
[ Info: Initializing simulation...
[ Info: Simulation time: 0 seconds, max(|uh|, |vh|, |h|): 6.28e-01, 0.00e+00, 1.00e+00 
[ Info:     ... simulation initialization complete (4.342 seconds)
[ Info: Executing initial time step...
[ Info:     ... initial time step complete (5.808 seconds).
[ Info: Simulation time: 990.000 ms, max(|uh|, |vh|, |h|): 3.45e-01, 5.25e-01, 1.00e+00 
[ Info: Simulation time: 1.990 seconds, max(|uh|, |vh|, |h|): 2.56e-01, 5.74e-01, 1.00e+00 
[ Info: Simulation time: 2.900 seconds, max(|uh|, |vh|, |h|): 6.10e-01, 1.50e-01, 1.00e+00 
[ Info: Simulation time: 3.810 seconds, max(|uh|, |vh|, |h|): 4.93e-01, 3.89e-01, 1.00e+00 
[ Info: Simulation time: 4.720 seconds, max(|uh|, |vh|, |h|): 4.78e-03, 6.28e-01, 1.00e+00 
[ Info: Simulation time: 5.630 seconds, max(|uh|, |vh|, |h|): 4.99e-01, 3.82e-01, 1.00e+00 
[ Info: Simulation time: 6.540 seconds, max(|uh|, |vh|, |h|): 6.08e-01, 1.60e-01, 1.00e+00 
[ Info: Simulation time: 7.450 seconds, max(|uh|, |vh|, |h|): 2.47e-01, 5.78e-01, 1.00e+00 
[ Info: Simulation time: 8.360 seconds, max(|uh|, |vh|, |h|): 3.05e-01, 5.50e-01, 1.00e+00 
[ Info: Simulation time: 9.270 seconds, max(|uh|, |vh|, |h|): 6.21e-01, 9.69e-02, 1.00e+00 
[ Info: Simulation time: 10.180 seconds, max(|uh|, |vh|, |h|): 4.57e-01, 4.31e-01, 1.00e+00 
[ Info: Simulation time: 11.090 seconds, max(|uh|, |vh|, |h|): 5.92e-02, 6.26e-01, 1.00e+00 
[ Info: Simulation time: 12.000 seconds, max(|uh|, |vh|, |h|): 5.30e-01, 3.37e-01, 1.00e+00 
[ Info: Simulation time: 12.900 seconds, max(|uh|, |vh|, |h|): 5.94e-01, 2.06e-01, 1.00e+00 
[ Info: Simulation time: 13.810 seconds, max(|uh|, |vh|, |h|): 2.02e-01, 5.95e-01, 1.00e+00 
[ Info: Simulation time: 14.720 seconds, max(|uh|, |vh|, |h|): 3.46e-01, 5.25e-01, 1.00e+00 
[ Info: Simulation time: 15.630 seconds, max(|uh|, |vh|, |h|): 6.26e-01, 4.89e-02, 1.00e+00 
[ Info: Simulation time: 16.590 seconds, max(|uh|, |vh|, |h|): 3.99e-01, 4.85e-01, 1.00e+00 
[ Info: Simulation time: 17.590 seconds, max(|uh|, |vh|, |h|): 1.92e-01, 5.98e-01, 1.00e+00 
[ Info: Simulation time: 18.590 seconds, max(|uh|, |vh|, |h|): 6.07e-01, 1.61e-01, 1.00e+00 
[ Info: Simulation time: 19.590 seconds, max(|uh|, |vh|, |h|): 4.64e-01, 4.24e-01, 1.00e+00 
[ Info: Simulation is stopping after running for 1.258 minutes.
[ Info: Simulation time 20 seconds equals or exceeds stop time 20 seconds.

````

## Perform Lagrangian filtering
Now we set up and run the offline Lagrangian filter on the output of the above simulation.
This could be performed in a different script (with appropriate import of Oceananigans.Units and CUDA if needed)

````julia
using OceananigansLagrangianFilter
````

### Set up the filter configuration

````julia
filter_config = OfflineFilterConfig(original_data_filename="SW_IO_with_tracer.jld2", # Where the original simulation output is
                                    output_filename = "SW_IO_with_tracer_filtered.jld2", # Where to save the filtered output
                                    var_names_to_filter = ("T",), # Which variables to filter
                                    velocity_names = ("u","v"), # Velocities to use for remapping
                                    architecture = CPU(), # CPU() or GPU(), if GPU() make sure you have CUDA.jl installed and imported
                                    Δt = 1e-2, # Time step of filtering simulation
                                    T_out=0.1, # How often to output filtered data
                                    N=2, # Order of Butterworth filter
                                    freq_c = 0.5, # Cut-off frequency of Butterworth filter
                                    compute_mean_velocities= true, # Whether to compute mean velocities
                                    output_netcdf = true, # Whether to output filtered data to a netcdf file in addition to .jld2
                                    delete_intermediate_files = true, # Delete the individual output of the forward and backward passes
                                    compute_Eulerian_filter = true) # Whether to compute the Eulerian filter for comparison
````

````
OfflineFilterConfig("SW_IO_with_tracer.jld2", ("T",), ("u", "v"), 0.0, 20.0, 20.0, CPU(), 0.1, (a1 = 0.17677669529663687, b1 = 0.1767766952966369, c1 = 0.35355339059327373, d1 = 0.3535533905932738, N_coeffs = 1), 0.01, InMemory{Int64}(1, 4), true, "forward_output.jld2", "backward_output.jld2", "SW_IO_with_tracer_filtered.jld2", 5, true, true, true, true, true, WENO{3, Float64, Float32}(order=5)
├── buffer_scheme: WENO{2, Float64, Float32}(order=3)
│   └── buffer_scheme: Centered(order=2)
└── advecting_velocity_scheme: Centered(order=4), 50×50×1 RectilinearGrid{Float64, Periodic, Periodic, Flat} on CPU with 3×3×0 halo
├── Periodic x ∈ [4.03717e-17, 6.28319) regularly spaced with Δx=0.125664
├── Periodic y ∈ [4.03717e-17, 6.28319) regularly spaced with Δy=0.125664
└── Flat z                              , "offline", "")
````

### Run the offline Lagrangian filter

````julia
run_offline_Lagrangian_filter(filter_config)
````

````
[ Info: Loaded data from SW_IO_with_tracer.jld2
[ Info: Created original variables: (:T,)
[ Info: Created filtered variables: (:T_C1, :xi_u_C1, :xi_v_C1, :T_S1, :xi_u_S1, :xi_v_S1)
[ Info: Created forcing for filtered variables
[ Info: Created model
[ Info: Initialised filtered variables
[ Info: Defined outputs
[ Info: Defined simulation
[ Info: Initializing simulation...
[ Info: Simulation time: 0 seconds
[ Info:     ... simulation initialization complete (1.545 minutes)
[ Info: Executing initial time step...
[ Info:     ... initial time step complete (19.100 seconds).
[ Info: Simulation time: 2 seconds
[ Info: Simulation time: 4 seconds
[ Info: Simulation time: 6 seconds
[ Info: Simulation time: 8 seconds
[ Info: Simulation time: 10 seconds
[ Info: Simulation time: 12 seconds
[ Info: Simulation time: 14 seconds
[ Info: Simulation time: 16 seconds
[ Info: Simulation time: 18 seconds
[ Info: Simulation is stopping after running for 5.971 minutes.
[ Info: Simulation time 20 seconds equals or exceeds stop time 20 seconds.
[ Info: Simulation time: 20 seconds
[ Info: Initializing simulation...
[ Info: Simulation time: 0 seconds
[ Info:     ... simulation initialization complete (112.529 ms)
[ Info: Executing initial time step...
[ Info:     ... initial time step complete (148.799 ms).
[ Info: Simulation time: 2 seconds
[ Info: Simulation time: 4 seconds
[ Info: Simulation time: 6 seconds
[ Info: Simulation time: 8 seconds
[ Info: Simulation time: 10 seconds
[ Info: Simulation time: 12 seconds
[ Info: Simulation time: 14 seconds
[ Info: Simulation time: 16 seconds
[ Info: Simulation time: 18 seconds
[ Info: Simulation is stopping after running for 4.049 minutes.
[ Info: Simulation time 20 seconds equals or exceeds stop time 20 seconds.
[ Info: Simulation time: 20 seconds
[ Info: Combined forward and backward contributions into SW_IO_with_tracer_filtered.jld2
[ Info: Wrote regridded data to new variables with _at_mean suffix in file SW_IO_with_tracer_filtered.jld2
[ Info: Computing Eulerian filter for variable T
[ Info: Computing Eulerian filter for variable u
[ Info: Computing Eulerian filter for variable v
[ Info: Wrote NetCDF file to SW_IO_with_tracer_filtered.nc

````

### Visualisation

````julia
using CairoMakie
````

Now we animate the results.

````julia
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
````

````
"IO_filtered_tracer_movie_offline.mp4"
````

![](IO_filtered_tracer_movie_offline.mp4)

We see that the Eulerian filter smudges the tracer field as it is advected by the
inertial oscillations. The Lagrangian means directly calculated by this method are
identical to the raw fields for the tracer and buoyancy shown, as they are conservative
fields. However, when we remap to a mean reference position, we see the value of the
Lagrangian filter in effectively removing the oscillations while preserving the tracer
structures.

````julia
# We remove these files to keep things tidy, keep them for analysis if desired
rm(filename_stem * ".jld2")
rm(filter_config.output_filename)
rm(filter_config.output_filename[1:end-5] * ".nc")
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


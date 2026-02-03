```@meta
EditURL = "../../../examples/offline_filter_lee_wave.jl"
```

# Steady lee waves with offline Lagrangian filtering

We set up a two-dimensional hydrostatic simulation of a steady flow over a
Gaussian bump, generating lee waves. We then apply an offline Lagrangian filter
to the output of the simulation to extract the mean flow, removing the wave
oscillations.

The simulation is initialised with a constant and uniform stratification and horizontal background velocity.

In this example, the filtering is performed offline after the simulation. We let the simulation run to steady
state, and run the filter on the steady data. This could also be performed online during the simulation, as
demonstrated in the `online_filter_geostrophic_adjustment` example.

## Run the simulation
### Install dependencies

````julia
using Oceananigans
using Oceananigans.Units
using NCDatasets
using Printf
using Oceananigans.Grids: xnode, znode
````

### Model parameters

Geometry

````julia
#Nx, Nz = 300, 400         # Higher res
Nx, Nz = 100, 100         # Lower res
H = 2kilometers           # Depth
L = 20kilometers          # Domain length
````

````
20000.0
````

Gaussian bump

````julia
h0 = 180meters            # Bump height
hill_width = 2kilometers  # Gaussian bump width
x0 = -10kilometers        # Bump location
hill(x) = h0 * exp(-(x-x0)^2 / 2hill_width^2)
bottom(x) = - H + hill(x)
````

````
bottom (generic function with 1 method)
````

Flow

````julia
f = 1e-4                  # Coriolis frequency [s⁻¹]
U = 0.2                   # Background flow speed [m/s]
Nsqr = 1e-6               # Buoyancy frequency squared [s⁻²]

κ = 0.1 # Diffusivity

filename_stem = "lee_wave";
````

### Define the grid

````julia
underlying_grid = RectilinearGrid(CPU(), size = (Nx, Nz), halo = (4, 4),
                                  x = (-L, L), z = (-H, 0),
                                  topology = (Periodic, Flat, Bounded))

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom))
````

````
100×1×100 ImmersedBoundaryGrid{Float64, Periodic, Flat, Bounded} on CPU with 4×0×4 halo:
├── immersed_boundary: PartialCellBottom(mean(zb)=-1979.4, min(zb)=-2000.0, max(zb)=-1824.0, ϵ=0.2)
├── underlying_grid: 100×1×100 RectilinearGrid{Float64, Periodic, Flat, Bounded} on CPU with 4×0×4 halo
├── Periodic x ∈ [-20000.0, 20000.0) regularly spaced with Δx=400.0
├── Flat y                           
└── Bounded  z ∈ [-2000.0, 0.0]      regularly spaced with Δz=20.0
````

### Closures

````julia
κ = 0.1
closure = ScalarDiffusivity(ν=κ, κ=κ)
````

````
ScalarDiffusivity{ExplicitTimeDiscretization}(ν=0.1, κ=0.1)
````

### Forcing

Steady geostrophic forcing through the v-momentum equation

````julia
A = f * U
@inline steady_forcing(x, z, t, p) = p.A
v_steady_forcing = Forcing(steady_forcing, parameters=(; A))
````

````
ContinuousForcing{@NamedTuple{A::Float64}}
├── func: steady_forcing (generic function with 1 method)
├── parameters: (A = 2.0e-5,)
└── field dependencies: ()
````

Sponge layer to absorb waves at periodic boundaries

````julia
t_restore = 4hour
sponge_x1 = -4*L/5
sponge_x2 = 4*L/5
sponge_dx = L/5

sponge_mask(x, p) = 1.0 + 0.5*(tanh((x + p.sponge_x1)/p.sponge_dx) - tanh((x + p.sponge_x2)/p.sponge_dx))
````

````
sponge_mask (generic function with 1 method)
````

HydrostaticFreeSurfaceModel has issues with continuous forcings on GPU, so we use discrete_form=true

````julia
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
````

````
DiscreteForcing{@NamedTuple{Nsqr::Float64, t_restore::Float64, sponge_x1::Float64, sponge_x2::Float64, sponge_dx::Float64}}
├── func: b_sponge_func (generic function with 1 method)
└── parameters: (Nsqr = 1.0e-6, t_restore = 14400.0, sponge_x1 = -16000.0, sponge_x2 = 16000.0, sponge_dx = 4000.0)
````

### Define the model

````julia
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
````

````
HydrostaticFreeSurfaceModel{CPU, ImmersedBoundaryGrid}(time = 0 seconds, iteration = 0)
├── grid: 100×1×100 ImmersedBoundaryGrid{Float64, Periodic, Flat, Bounded} on CPU with 4×0×4 halo
├── timestepper: QuasiAdamsBashforth2TimeStepper
├── tracers: b
├── closure: ScalarDiffusivity{ExplicitTimeDiscretization}(ν=0.1, κ=(b=0.1,))
├── buoyancy: BuoyancyTracer with ĝ = NegativeZDirection()
├── free surface: SplitExplicitFreeSurface with gravitational acceleration 9.80665 m s⁻²
│   └── substepping: FixedTimeStepSize(1.999 seconds)
├── advection scheme: 
│   ├── momentum: WENO{3, Float64, Float32}(order=5)
│   └── b: WENO{3, Float64, Float32}(order=5)
├── vertical_coordinate: ZCoordinate
└── coriolis: FPlane{Float64}
````

### Initialise the buoyancy

````julia
bᵢ(x, z) = Nsqr * z
set!(model, u=U, b=bᵢ)
````

### Define the simulation

````julia
simulation = Simulation(model, Δt=20seconds, stop_time=15days)
````

````
Simulation of HydrostaticFreeSurfaceModel{CPU, ImmersedBoundaryGrid}(time = 0 seconds, iteration = 0)
├── Next time step: 20 seconds
├── run_wall_time: 0 seconds
├── run_wall_time / iteration: NaN days
├── stop_time: 15 days
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

Set an adaptive timestep

````julia
conjure_time_step_wizard!(simulation, IterationInterval(500), cfl=0.2, max_Δt=2minutes)
````

Add a progress callback

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

add_callback!(simulation, print_progress, IterationInterval(500))
````

### Set up the output

````julia
u, v, w = model.velocities
b = model.tracers.b
````

````
100×1×100 Field{Center, Center, Center} on ImmersedBoundaryGrid on CPU
├── grid: 100×1×100 ImmersedBoundaryGrid{Float64, Periodic, Flat, Bounded} on CPU with 4×0×4 halo
├── boundary conditions: FieldBoundaryConditions
│   └── west: Periodic, east: Periodic, south: Nothing, north: Nothing, bottom: ZeroFlux, top: ZeroFlux, immersed: ZeroFlux
└── data: 108×1×108 OffsetArray(::Array{Float64, 3}, -3:104, 1:1, -3:104) with eltype Float64 with indices -3:104×1:1×-3:104
    └── max=-1.0e-5, min=-0.00199, mean=-0.000990905
````

Output a jld2 file for Lagrangian filtering

````julia
simulation.output_writers[:jld2fields] = JLD2Writer(
    model, (; u, w, b), filename = filename_stem * ".jld2", schedule=TimeInterval(1hour), overwrite_existing=true)
````

````
JLD2Writer scheduled on TimeInterval(1 hour):
├── filepath: lee_wave.jld2
├── 3 outputs: (u, w, b)
├── array_type: Array{Float32}
├── including: [:grid, :coriolis, :buoyancy, :closure]
├── file_splitting: NoFileSplitting
└── file size: 0 bytes (file not yet created)
````

### Run simulation

````julia
@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)
````

````
[ Info: Running the simulation...
[ Info: Initializing simulation...
[00.00%] i: 0, t: 0 seconds, wall time: 7.316 seconds, max(u): (2.000e-01, 0.000e+00, 1.203e-02) m/s, next Δt: 22 seconds
[ Info:     ... simulation initialization complete (5.625 seconds)
[ Info: Executing initial time step...
[ Info:     ... initial time step complete (23.988 seconds).
[00.85%] i: 500, t: 3.049 hours, wall time: 30.393 seconds, max(u): (3.600e-01, 6.641e-02, 1.477e-02) m/s, next Δt: 24.200 seconds
[01.78%] i: 1000, t: 6.403 hours, wall time: 3.029 seconds, max(u): (3.909e-01, 1.453e-01, 1.707e-02) m/s, next Δt: 26.620 seconds
[02.80%] i: 1500, t: 10.081 hours, wall time: 2.791 seconds, max(u): (4.185e-01, 1.604e-01, 1.845e-02) m/s, next Δt: 29.282 seconds
[03.93%] i: 2000, t: 14.146 hours, wall time: 2.855 seconds, max(u): (4.259e-01, 1.522e-01, 1.921e-02) m/s, next Δt: 32.210 seconds
[05.17%] i: 2500, t: 18.608 hours, wall time: 2.866 seconds, max(u): (4.223e-01, 1.477e-01, 1.990e-02) m/s, next Δt: 35.431 seconds
[06.53%] i: 3000, t: 23.512 hours, wall time: 2.903 seconds, max(u): (4.247e-01, 1.363e-01, 1.996e-02) m/s, next Δt: 38.974 seconds
[08.02%] i: 3500, t: 1.204 days, wall time: 2.920 seconds, max(u): (4.414e-01, 1.321e-01, 2.056e-02) m/s, next Δt: 42.872 seconds
[09.68%] i: 4000, t: 1.451 days, wall time: 2.702 seconds, max(u): (4.477e-01, 1.331e-01, 2.059e-02) m/s, next Δt: 47.159 seconds
[11.48%] i: 4500, t: 1.722 days, wall time: 3.294 seconds, max(u): (4.404e-01, 1.296e-01, 2.027e-02) m/s, next Δt: 51.875 seconds
[13.47%] i: 5000, t: 2.020 days, wall time: 2.776 seconds, max(u): (4.362e-01, 1.259e-01, 1.993e-02) m/s, next Δt: 57.062 seconds
[15.63%] i: 5500, t: 2.345 days, wall time: 3.060 seconds, max(u): (4.358e-01, 1.291e-01, 2.031e-02) m/s, next Δt: 1.046 minutes
[18.03%] i: 6000, t: 2.705 days, wall time: 2.882 seconds, max(u): (4.446e-01, 1.334e-01, 2.027e-02) m/s, next Δt: 1.151 minutes
[20.66%] i: 6500, t: 3.099 days, wall time: 3.550 seconds, max(u): (4.293e-01, 1.287e-01, 1.969e-02) m/s, next Δt: 1.266 minutes
[23.55%] i: 7000, t: 3.533 days, wall time: 3.235 seconds, max(u): (4.423e-01, 1.287e-01, 2.051e-02) m/s, next Δt: 1.392 minutes
[26.71%] i: 7500, t: 4.006 days, wall time: 3.387 seconds, max(u): (4.355e-01, 1.287e-01, 1.984e-02) m/s, next Δt: 1.532 minutes
[30.18%] i: 8000, t: 4.528 days, wall time: 3.159 seconds, max(u): (4.374e-01, 1.286e-01, 2.026e-02) m/s, next Δt: 1.656 minutes
[33.93%] i: 8500, t: 5.090 days, wall time: 3.523 seconds, max(u): (4.355e-01, 1.288e-01, 1.982e-02) m/s, next Δt: 1.670 minutes
[37.79%] i: 9000, t: 5.669 days, wall time: 3.633 seconds, max(u): (4.408e-01, 1.297e-01, 2.041e-02) m/s, next Δt: 1.644 minutes
[41.55%] i: 9500, t: 6.232 days, wall time: 3.828 seconds, max(u): (4.332e-01, 1.296e-01, 1.986e-02) m/s, next Δt: 1.675 minutes
[45.40%] i: 10000, t: 6.810 days, wall time: 3.234 seconds, max(u): (4.416e-01, 1.301e-01, 2.032e-02) m/s, next Δt: 1.643 minutes
[49.16%] i: 10500, t: 7.373 days, wall time: 3.104 seconds, max(u): (4.339e-01, 1.297e-01, 1.999e-02) m/s, next Δt: 1.671 minutes
[53.01%] i: 11000, t: 7.951 days, wall time: 3.248 seconds, max(u): (4.402e-01, 1.304e-01, 2.023e-02) m/s, next Δt: 1.649 minutes
[56.77%] i: 11500, t: 8.515 days, wall time: 2.817 seconds, max(u): (4.353e-01, 1.296e-01, 2.007e-02) m/s, next Δt: 1.666 minutes
[60.52%] i: 12000, t: 9.079 days, wall time: 2.835 seconds, max(u): (4.392e-01, 1.294e-01, 2.021e-02) m/s, next Δt: 1.653 minutes
[64.27%] i: 12500, t: 9.641 days, wall time: 2.851 seconds, max(u): (4.358e-01, 1.285e-01, 2.006e-02) m/s, next Δt: 1.665 minutes
[68.03%] i: 13000, t: 10.205 days, wall time: 2.811 seconds, max(u): (4.386e-01, 1.289e-01, 2.015e-02) m/s, next Δt: 1.656 minutes
[71.78%] i: 13500, t: 10.767 days, wall time: 2.847 seconds, max(u): (4.360e-01, 1.285e-01, 2.007e-02) m/s, next Δt: 1.664 minutes
[75.54%] i: 14000, t: 11.331 days, wall time: 2.812 seconds, max(u): (4.381e-01, 1.287e-01, 2.011e-02) m/s, next Δt: 1.657 minutes
[79.29%] i: 14500, t: 11.893 days, wall time: 2.834 seconds, max(u): (4.365e-01, 1.285e-01, 2.011e-02) m/s, next Δt: 1.662 minutes
[83.05%] i: 15000, t: 12.457 days, wall time: 2.813 seconds, max(u): (4.375e-01, 1.287e-01, 2.008e-02) m/s, next Δt: 1.660 minutes
[86.80%] i: 15500, t: 13.020 days, wall time: 2.836 seconds, max(u): (4.371e-01, 1.285e-01, 2.012e-02) m/s, next Δt: 1.660 minutes
[90.55%] i: 16000, t: 13.583 days, wall time: 2.819 seconds, max(u): (4.371e-01, 1.287e-01, 2.008e-02) m/s, next Δt: 1.661 minutes
[94.31%] i: 16500, t: 14.146 days, wall time: 2.835 seconds, max(u): (4.373e-01, 1.286e-01, 2.012e-02) m/s, next Δt: 1.660 minutes
[98.06%] i: 17000, t: 14.708 days, wall time: 2.808 seconds, max(u): (4.370e-01, 1.286e-01, 2.008e-02) m/s, next Δt: 1.661 minutes
[ Info: Simulation is stopping after running for 2.242 minutes.
[ Info: Simulation time 15 days equals or exceeds stop time 15 days.
[ Info: Simulation completed in 2.242 minutes

````

## Perform Lagrangian filtering
Now we set up and run the offline Lagrangian filter on the output of the above simulation.
This could be performed in a different script (with appropriate import of Oceananigans.Units and CUDA if needed)

````julia
using OceananigansLagrangianFilter
````

### Set up the filter configuration

````julia
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
````

````
OfflineFilterConfig("lee_wave.jld2", ("b",), ("u", "w"), 432000.0, 1.296e6, 864000.0, CPU(), 3600.0, (a1 = 4.783542904563622e-6, b1 = 1.1548494156391084e-5, c1 = 1.913417161825449e-5, d1 = 4.619397662556434e-5, a2 = 1.1548494156391084e-5, b2 = 4.783542904563623e-6, c2 = 4.619397662556434e-5, d2 = 1.9134171618254493e-5, N_coeffs = 2), 60.0, InMemory{Int64}(1, 4), true, "forward_output.jld2", "backward_output.jld2", "lee_wave_offline_filtered.jld2", 5, true, true, true, true, true, WENO{3, Float64, Float32}(order=5)
├── buffer_scheme: WENO{2, Float64, Float32}(order=3)
│   └── buffer_scheme: Centered(order=2)
└── advecting_velocity_scheme: Centered(order=4), 100×1×100 ImmersedBoundaryGrid{Float64, Periodic, Flat, Bounded} on CPU with 4×0×4 halo:
├── immersed_boundary: PartialCellBottom(mean(zb)=-1979.4, min(zb)=-2000.0, max(zb)=-1824.0, ϵ=0.2)
├── underlying_grid: 100×1×100 RectilinearGrid{Float64, Periodic, Flat, Bounded} on CPU with 4×0×4 halo
├── Periodic x ∈ [-20000.0, 20000.0) regularly spaced with Δx=400.0
├── Flat y                           
└── Bounded  z ∈ [-2000.0, 0.0]      regularly spaced with Δz=20.0, "offline", "")
````

### Run the offline Lagrangian filter

````julia
run_offline_Lagrangian_filter(filter_config)
````

````
[ Info: Loaded data from lee_wave.jld2
[ Info: Created original variables: (:b,)
[ Info: Created filtered variables: (:b_C1, :b_C2, :xi_u_C1, :xi_u_C2, :xi_w_C1, :xi_w_C2, :b_S1, :b_S2, :xi_u_S1, :xi_u_S2, :xi_w_S1, :xi_w_S2)
[ Info: Created forcing for filtered variables
[ Info: Created model
[ Info: Initialised filtered variables
[ Info: Defined outputs
[ Info: Defined simulation
[ Info: Initializing simulation...
[ Info: Simulation time: 0 seconds
[ Info:     ... simulation initialization complete (24.134 minutes)
[ Info: Executing initial time step...
[ Info:     ... initial time step complete (1.253 minutes).
[ Info: Simulation time: 1 day
[ Info: Simulation time: 2 days
[ Info: Simulation time: 3 days
[ Info: Simulation time: 4 days
[ Info: Simulation time: 5 days
[ Info: Simulation time: 6 days
[ Info: Simulation time: 7 days
[ Info: Simulation time: 8 days
[ Info: Simulation time: 9 days
[ Info: Simulation is stopping after running for 6.889 hours.
[ Info: Simulation time 10 days equals or exceeds stop time 10 days.
[ Info: Simulation time: 10 days
[ Info: Initializing simulation...
[ Info: Simulation time: 0 seconds
[ Info:     ... simulation initialization complete (649.002 ms)
[ Info: Executing initial time step...
[ Info:     ... initial time step complete (2.226 seconds).
[ Info: Simulation time: 1 day
[ Info: Simulation time: 2 days
[ Info: Simulation time: 3 days
[ Info: Simulation time: 4 days
[ Info: Simulation time: 5 days
[ Info: Simulation time: 6 days
[ Info: Simulation time: 7 days
[ Info: Simulation time: 8 days
[ Info: Simulation time: 9 days
[ Info: Simulation is stopping after running for 6.491 hours.
[ Info: Simulation time 10 days equals or exceeds stop time 10 days.
[ Info: Simulation time: 10 days
[ Info: Combined forward and backward contributions into lee_wave_offline_filtered.jld2
┌ Warning: Grid is an ImmersedBoundaryGrid, regridding to mean position may not be sensible. 
│             Areas below bottom height will be masked with NaNs, assuming bottom height is a function of x and/or y.
└ @ OceananigansLagrangianFilter.Utils ~/Documents/Projects/OceananigansLagrangianFilter/src/Utils/post_processing_utils.jl:357
[ Info: Wrote regridded data to new variables with _at_mean suffix in file lee_wave_offline_filtered.jld2
[ Info: Computing Eulerian filter for variable b
[ Info: Computing Eulerian filter for variable u
[ Info: Computing Eulerian filter for variable w
[ Info: Wrote NetCDF file to lee_wave_offline_filtered.nc

````

### Visualisation

````julia
using CairoMakie
````

Now we animate the filterieng results for the horizontal velocity u and buoyancy b:

````julia
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
````

````
Makie.Colorbar()
````

Plot topography

````julia
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
````

````
"lee_wave_filtered_u_movie_offline.mp4"
````

![](lee_wave_filtered_u_movie_offline.mp4)

The filter has been run after the initial simulation spin-up (time shown is time from day 5 of the
original simulation, so the raw fields are already in steady state). The filtered fields also need
time to spin-up (corresponding to the characteristic timescale of the cut-off filter, here ~1.5 days),
so the initial and final frames of the animation show transient behaviour in the filtered fields.

The Eulerian filtered output looks very similar to the raw output, since the lee waves are steady.
They're steady because they've been Doppler-shifted by the mean flow. When we use the Lagrangian
filter instead, we see that the lee waves are removed as they are high frequency in the frame of the
flow. Note that the colour range shown is an order of magnitude smaller for the Lagrangian filtered
fields, allowing us to see the wave impact on the mean flow.

The contours of buoyancy demonstrate the difference between the displaced Lagrangian filtered field
(bottom left) and the Lagrangian filtered field remapped to the mean position (bottom right). This
removes the wave-like displacements from the isopycnals. The interpolation stage to calculate this
field is imperfect near the immersed boundary (it's not guaranteed that every spatial location has
a corresponding Lagrangian mean defined at the mean position), so it has been masked.

We remove these files to keep things tidy, keep them for analysis if desired

````julia
rm(filename_stem * ".jld2")
rm(filter_config.output_filename)
rm(filter_config.output_filename[1:end-5] * ".nc")
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


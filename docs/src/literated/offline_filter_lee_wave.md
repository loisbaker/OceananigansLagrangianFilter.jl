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

κv = 0.1 # Vertical diffusivity
κBh = L^2 / pi^2/ Nx^2 * 100 # Horizontal biharmonic diffusivity

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
100×1×100 ImmersedBoundaryGrid{Float64, Oceananigans.Grids.Periodic, Oceananigans.Grids.Flat, Oceananigans.Grids.Bounded} on CPU with 4×0×4 halo:
├── immersed_boundary: PartialCellBottom(mean(zb)=-1979.4, min(zb)=-2000.0, max(zb)=-1824.0, ϵ=0.2)
├── underlying_grid: 100×1×100 RectilinearGrid{Float64, Periodic, Flat, Bounded} on CPU with 4×0×4 halo
├── Periodic x ∈ [-20000.0, 20000.0) regularly spaced with Δx=400.0
├── Flat y                           
└── Bounded  z ∈ [-2000.0, 0.0]      regularly spaced with Δz=20.0
````

### Closures

````julia
horizontal_closure = HorizontalScalarBiharmonicDiffusivity(ν=κBh, κ=κBh)
vertical_closure = VerticalScalarDiffusivity(ν=κv, κ=κv)
closure = (horizontal_closure, vertical_closure)
````

````
(ScalarBiharmonicDiffusivity{HorizontalFormulation}(ν=4152.85, κ=4152.85), VerticalScalarDiffusivity{ExplicitTimeDiscretization}(ν=0.1, κ=0.1))
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
                                    forcing = (; v = (v_steady_forcing, v_sponge_forcing),
                                                u = u_sponge_forcing, b = b_sponge_forcing),
                                    closure=closure,
                                    )
````

````
HydrostaticFreeSurfaceModel{CPU, ImmersedBoundaryGrid}(time = 0 seconds, iteration = 0)
├── grid: 100×1×100 ImmersedBoundaryGrid{Float64, Oceananigans.Grids.Periodic, Oceananigans.Grids.Flat, Oceananigans.Grids.Bounded} on CPU with 4×0×4 halo
├── timestepper: QuasiAdamsBashforth2TimeStepper
├── tracers: b
├── closure: Tuple with 2 closures:
│   ├── ScalarBiharmonicDiffusivity{HorizontalFormulation}(ν=4152.85, κ=(b=4152.85,))
│   └── VerticalScalarDiffusivity{ExplicitTimeDiscretization}(ν=0.1, κ=(b=0.1,))
├── buoyancy: Oceananigans.BuoyancyFormulations.BuoyancyTracer with ĝ = NegativeZDirection()
├── free surface: Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces.SplitExplicitFreeSurface with gravitational acceleration 9.80665 m s⁻²
│   └── substepping: FixedTimeStepSize(1.999 seconds)
├── advection scheme: 
│   ├── momentum: VectorInvariant
│   └── b: Centered(order=2)
├── vertical_coordinate: Oceananigans.Models.HydrostaticFreeSurfaceModels.ZCoordinate
└── coriolis: Oceananigans.Coriolis.FPlane{Oceananigans.Advection.EnstrophyConserving{Float64}, Float64}
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
100×1×100 Field{Oceananigans.Grids.Center, Oceananigans.Grids.Center, Oceananigans.Grids.Center} on Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid on CPU
├── grid: 100×1×100 ImmersedBoundaryGrid{Float64, Oceananigans.Grids.Periodic, Oceananigans.Grids.Flat, Oceananigans.Grids.Bounded} on CPU with 4×0×4 halo
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
├── including: [:coriolis, :buoyancy, :closure]
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
[00.00%] i: 0, t: 0 seconds, wall time: 7.916 seconds, max(u): (2.000e-01, 0.000e+00, 1.203e-02) m/s, next Δt: 22 seconds
[ Info:     ... simulation initialization complete (9.158 seconds)
[ Info: Executing initial time step...
[ Info:     ... initial time step complete (30.090 seconds).
[00.85%] i: 500, t: 3.049 hours, wall time: 38.104 seconds, max(u): (3.733e-01, 7.389e-02, 1.929e-02) m/s, next Δt: 24.200 seconds
[01.78%] i: 1000, t: 6.403 hours, wall time: 2.421 seconds, max(u): (4.306e-01, 1.541e-01, 1.885e-02) m/s, next Δt: 26.620 seconds
[02.80%] i: 1500, t: 10.081 hours, wall time: 2.129 seconds, max(u): (4.696e-01, 1.716e-01, 2.007e-02) m/s, next Δt: 29.282 seconds
[03.93%] i: 2000, t: 14.146 hours, wall time: 2.259 seconds, max(u): (4.843e-01, 1.758e-01, 2.162e-02) m/s, next Δt: 32.210 seconds
[05.17%] i: 2500, t: 18.608 hours, wall time: 2.210 seconds, max(u): (4.880e-01, 1.626e-01, 2.103e-02) m/s, next Δt: 35.431 seconds
[06.53%] i: 3000, t: 23.512 hours, wall time: 2.302 seconds, max(u): (4.919e-01, 1.674e-01, 2.090e-02) m/s, next Δt: 38.974 seconds
[08.02%] i: 3500, t: 1.204 days, wall time: 2.263 seconds, max(u): (5.099e-01, 1.763e-01, 2.151e-02) m/s, next Δt: 42.872 seconds
[09.68%] i: 4000, t: 1.451 days, wall time: 2.074 seconds, max(u): (5.114e-01, 1.782e-01, 2.152e-02) m/s, next Δt: 47.159 seconds
[11.48%] i: 4500, t: 1.722 days, wall time: 2.822 seconds, max(u): (5.038e-01, 1.748e-01, 2.107e-02) m/s, next Δt: 51.875 seconds
[13.47%] i: 5000, t: 2.020 days, wall time: 2.193 seconds, max(u): (4.963e-01, 1.733e-01, 2.075e-02) m/s, next Δt: 57.062 seconds
[15.63%] i: 5500, t: 2.345 days, wall time: 2.511 seconds, max(u): (5.005e-01, 1.719e-01, 2.094e-02) m/s, next Δt: 1.046 minutes
[18.03%] i: 6000, t: 2.705 days, wall time: 2.322 seconds, max(u): (5.009e-01, 1.731e-01, 2.108e-02) m/s, next Δt: 1.151 minutes
[20.66%] i: 6500, t: 3.099 days, wall time: 3.302 seconds, max(u): (4.899e-01, 1.703e-01, 2.099e-02) m/s, next Δt: 1.266 minutes
[23.55%] i: 7000, t: 3.533 days, wall time: 2.880 seconds, max(u): (5.061e-01, 1.737e-01, 2.145e-02) m/s, next Δt: 1.392 minutes
[26.71%] i: 7500, t: 4.006 days, wall time: 3.052 seconds, max(u): (5.007e-01, 1.715e-01, 2.085e-02) m/s, next Δt: 1.453 minutes
[30.01%] i: 8000, t: 4.502 days, wall time: 2.691 seconds, max(u): (5.028e-01, 1.729e-01, 2.081e-02) m/s, next Δt: 1.451 minutes
[33.32%] i: 8500, t: 4.999 days, wall time: 1.907 seconds, max(u): (4.986e-01, 1.717e-01, 2.087e-02) m/s, next Δt: 1.456 minutes
[36.63%] i: 9000, t: 5.495 days, wall time: 1.914 seconds, max(u): (4.997e-01, 1.713e-01, 2.090e-02) m/s, next Δt: 1.453 minutes
[39.94%] i: 9500, t: 5.991 days, wall time: 2.206 seconds, max(u): (5.012e-01, 1.723e-01, 2.089e-02) m/s, next Δt: 1.451 minutes
[43.24%] i: 10000, t: 6.487 days, wall time: 1.896 seconds, max(u): (4.995e-01, 1.717e-01, 2.081e-02) m/s, next Δt: 1.456 minutes
[46.55%] i: 10500, t: 6.983 days, wall time: 2.223 seconds, max(u): (5.009e-01, 1.724e-01, 2.091e-02) m/s, next Δt: 1.451 minutes
[49.86%] i: 11000, t: 7.478 days, wall time: 2.270 seconds, max(u): (4.986e-01, 1.715e-01, 2.083e-02) m/s, next Δt: 1.457 minutes
[53.16%] i: 11500, t: 7.975 days, wall time: 1.904 seconds, max(u): (5.011e-01, 1.723e-01, 2.091e-02) m/s, next Δt: 1.451 minutes
[56.47%] i: 12000, t: 8.470 days, wall time: 1.903 seconds, max(u): (4.990e-01, 1.717e-01, 2.082e-02) m/s, next Δt: 1.457 minutes
[59.78%] i: 12500, t: 8.966 days, wall time: 1.903 seconds, max(u): (5.008e-01, 1.722e-01, 2.089e-02) m/s, next Δt: 1.451 minutes
[63.08%] i: 13000, t: 9.462 days, wall time: 1.903 seconds, max(u): (4.991e-01, 1.717e-01, 2.084e-02) m/s, next Δt: 1.456 minutes
[66.39%] i: 13500, t: 9.958 days, wall time: 1.878 seconds, max(u): (5.005e-01, 1.721e-01, 2.088e-02) m/s, next Δt: 1.452 minutes
[69.70%] i: 14000, t: 10.455 days, wall time: 1.886 seconds, max(u): (4.996e-01, 1.719e-01, 2.084e-02) m/s, next Δt: 1.455 minutes
[73.01%] i: 14500, t: 10.951 days, wall time: 1.906 seconds, max(u): (5.001e-01, 1.720e-01, 2.087e-02) m/s, next Δt: 1.453 minutes
[76.31%] i: 15000, t: 11.447 days, wall time: 1.913 seconds, max(u): (4.998e-01, 1.719e-01, 2.086e-02) m/s, next Δt: 1.454 minutes
[79.62%] i: 15500, t: 11.943 days, wall time: 1.937 seconds, max(u): (5.000e-01, 1.719e-01, 2.086e-02) m/s, next Δt: 1.454 minutes
[82.93%] i: 16000, t: 12.439 days, wall time: 1.927 seconds, max(u): (5.000e-01, 1.720e-01, 2.086e-02) m/s, next Δt: 1.454 minutes
[86.23%] i: 16500, t: 12.935 days, wall time: 1.916 seconds, max(u): (4.998e-01, 1.719e-01, 2.085e-02) m/s, next Δt: 1.454 minutes
[89.54%] i: 17000, t: 13.431 days, wall time: 1.912 seconds, max(u): (5.001e-01, 1.720e-01, 2.087e-02) m/s, next Δt: 1.453 minutes
[92.85%] i: 17500, t: 13.927 days, wall time: 1.945 seconds, max(u): (4.998e-01, 1.719e-01, 2.086e-02) m/s, next Δt: 1.454 minutes
[96.15%] i: 18000, t: 14.423 days, wall time: 1.952 seconds, max(u): (5.001e-01, 1.720e-01, 2.087e-02) m/s, next Δt: 1.453 minutes
[99.46%] i: 18500, t: 14.919 days, wall time: 1.952 seconds, max(u): (4.998e-01, 1.719e-01, 2.086e-02) m/s, next Δt: 1.454 minutes
[ Info: Simulation is stopping after running for 1.999 minutes.
[ Info: Simulation time 15 days equals or exceeds stop time 15 days.
[ Info: Simulation completed in 1.999 minutes

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
OceananigansLagrangianFilter.OfflineLagrangianFilter.OfflineFilterConfig("lee_wave.jld2", ("b",), ("u", "w"), 432000.0, 1.296e6, 864000.0, Oceananigans.Architectures.CPU(), 3600.0, (a1 = 4.783542904563622e-6, b1 = 1.1548494156391084e-5, c1 = 1.913417161825449e-5, d1 = 4.619397662556434e-5, a2 = 1.1548494156391084e-5, b2 = 4.783542904563623e-6, c2 = 4.619397662556434e-5, d2 = 1.9134171618254493e-5, N_coeffs = 2), 60.0, Oceananigans.OutputReaders.InMemory{Int64}(1, 4), true, "forward_output.jld2", "backward_output.jld2", "lee_wave_offline_filtered.jld2", 5, true, true, true, true, true, WENO{3, Float64, Nothing}(order=5)
├── buffer_scheme: WENO{2, Float64, Nothing}(order=3)
│   └── buffer_scheme: Centered(order=2)
└── advecting_velocity_scheme: Centered(order=4), 100×1×100 ImmersedBoundaryGrid{Float64, Oceananigans.Grids.Periodic, Oceananigans.Grids.Flat, Oceananigans.Grids.Bounded} on CPU with 4×0×4 halo:
├── immersed_boundary: PartialCellBottom(mean(zb)=-1979.4, min(zb)=-2000.0, max(zb)=-1824.0, ϵ=0.2)
├── underlying_grid: 100×1×100 RectilinearGrid{Float64, Periodic, Flat, Bounded} on CPU with 4×0×4 halo
├── Periodic x ∈ [-20000.0, 20000.0) regularly spaced with Δx=400.0
├── Flat y                           
└── Bounded  z ∈ [-2000.0, 0.0]      regularly spaced with Δz=20.0, "", false, nothing, nothing, nothing)
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
[ Info:     ... simulation initialization complete (24.933 minutes)
[ Info: Executing initial time step...
[ Info:     ... initial time step complete (32.720 seconds).
[ Info: Simulation time: 1 day
[ Info: Simulation time: 2 days
[ Info: Simulation time: 3 days
[ Info: Simulation time: 4 days
[ Info: Simulation time: 5 days
[ Info: Simulation time: 6 days
[ Info: Simulation time: 7 days
[ Info: Simulation time: 8 days
[ Info: Simulation time: 9 days
[ Info: Simulation is stopping after running for 6.578 hours.
[ Info: Simulation time 10 days equals or exceeds stop time 10 days.
[ Info: Simulation time: 10 days
[ Info: Initializing simulation...
[ Info: Simulation time: 0 seconds
[ Info:     ... simulation initialization complete (695.246 ms)
[ Info: Executing initial time step...
[ Info:     ... initial time step complete (2.206 seconds).
[ Info: Simulation time: 1 day
[ Info: Simulation time: 2 days
[ Info: Simulation time: 3 days
[ Info: Simulation time: 4 days
[ Info: Simulation time: 5 days
[ Info: Simulation time: 6 days
[ Info: Simulation time: 7 days
[ Info: Simulation time: 8 days
[ Info: Simulation time: 9 days
[ Info: Simulation is stopping after running for 6.149 hours.
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


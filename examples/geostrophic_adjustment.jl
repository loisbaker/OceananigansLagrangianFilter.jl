using Oceananigans 
using Oceananigans.Units
using Oceananigans.TurbulenceClosures
using CairoMakie 
using NCDatasets
using Printf

# We set up a geostrophic adjustment problem similar to Blumen (2000), JPO
# in a domain that is horizontally periodic.

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

# No real need for a closure here, but we include one for demo with the online filter
horizontal_closure = HorizontalScalarDiffusivity(ν=κh, κ=κh )
vertical_closure = VerticalScalarDiffusivity(ν=κv , κ=κv )
closure = (horizontal_closure, vertical_closure)

# Define tracers
tracers = (:b,:T)

# Define the model
model =  NonhydrostaticModel(; grid,
                coriolis = FPlane(f = f),
                buoyancy = BuoyancyTracer(),
                tracers = tracers,
                advection = WENO(),
                closure = closure)

# Initialise the buoyancy and tracer (velocities start at rest by default)
bᵢ(x, z) = Δb*sin(2*pi/L_front * x)
Tᵢ(x, z) = exp(-(x/(L_front/50)).^2)
set!(model, b= bᵢ, T= Tᵢ) 

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
u, v, w = model.velocities
b = model.tracers.b
T = model.tracers.T
# Fields to be filtered must be specified at cell centers, so we can interpolate before 
# output if necessary. 
wc = Field(@at (Center, Center, Center) model.velocities.w) 

# For Lagrangian filtering
simulation.output_writers[:jld2fields] = JLD2Writer(
    model, (; b, u, v, w, wc, T), filename = filename_stem * ".jld2", schedule=TimeInterval(1hour), overwrite_existing=true)


# NetCDF can be useful too for visualisation
rm(filename_stem * ".nc",force=true)
simulation.output_writers[:ncfields] = NetCDFWriter(
    model, (; b, u, v, w, wc,T), filename = filename_stem * ".nc", schedule=TimeInterval(1hour), overwrite_existing=true)
    
@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)
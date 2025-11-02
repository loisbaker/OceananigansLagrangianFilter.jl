# # Geostrophic adjustment in 3D with offline Lagrangian filtering

# We set up a geostrophic adjustment problem similar to Blumen (2000), JPO
# in a domain that is horizontally periodic. We include a few gridpoints in the along-front
# direction (y) to demonstrate 3D capability.

# An initially unbalanced two-dimensional
# front oscillates with the inertial frequency around a state of geostrophic balance,
# and we illustrate that we can remove the oscillations to find the mean state.
# Credit to Tom Cummings for work on this example. 

# In this example, the filtering is performed offline after the simulation.

using Oceananigans 
using Oceananigans.Units
using NCDatasets
using Printf

# We set up a geostrophic adjustment problem similar to Blumen (2000), JPO
# in a domain that is horizontally periodic. Credit to Tom Cummings for work on this example. 

# Model parameters
Nx = 120
Nz = 40
Ny = 5 # Just a few points in y to make it 3D, but not too computationally expensive
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

filename_stem = "geostrophic_adjustment_3D"

grid = RectilinearGrid(CPU(),size = (Nx, Ny, Nz), 
                       x = (-L_front/2, L_front/2),
                       z = (-H, 0),
                       y = (0, L_front/Nx*Ny),
                       topology = (Periodic, Periodic, Bounded))

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
bᵢ(x, y, z) = Δb*sin(2*pi/L_front * x)
Tᵢ(x, y, z) = exp(-x.^2/(L_front/50).^2)
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

# For Lagrangian filtering
simulation.output_writers[:jld2fields] = JLD2Writer(
    model, (; b, u, v, w, T), filename = filename_stem * ".jld2", schedule=TimeInterval(1hour), overwrite_existing=true)


# NetCDF can be useful too for visualisation
rm(filename_stem * ".nc",force=true)
simulation.output_writers[:ncfields] = NetCDFWriter(
    model, (; b, u, v, w, T), filename = filename_stem * ".nc", schedule=TimeInterval(1hour), overwrite_existing=true)
    
@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)

# Now we set up and run the offline Lagrangian filter on the output of the above simulation.
# This could be performed in a different script (with appropriate import of Oceananigans.Units and CUDA if needed)

using OceananigansLagrangianFilter


# Set up the filter configuration
filter_config = OfflineFilterConfig(original_data_filename="geostrophic_adjustment_3D.jld2", # Where the original simulation output is
                                    output_filename = "geostrophic_adjustment_offline_filtered_3D.jld2", # Where to save the filtered output
                                    var_names_to_filter = ("T", "b"), # Which variables to filter
                                    velocity_names = ("u","v","w"), # Velocities to use for Lagrangian filtering
                                    architecture = CPU(), # CPU() or GPU(), if GPU() make sure you have CUDA.jl installed and imported
                                    Δt = 20minutes, # Time step of filtering simulation
                                    T_out = 1hour, # How often to output filtered data
                                    N = 2, # Order of Butterworth filter
                                    freq_c = 1e-4/2, # Cut-off frequency of Butterworth filter
                                    compute_mean_velocities= true, # Whether to compute mean velocities 
                                    output_netcdf = true, # Whether to output filtered data to a netcdf file in addition to .jld2
                                    delete_intermediate_files = true, # Delete the individual output of the forward and backward passes
                                    compute_Eulerian_filter = true) # Whether to compute the Eulerian filter for comparison


# Run the offline Lagrangian filter
run_offline_Lagrangian_filter(filter_config)


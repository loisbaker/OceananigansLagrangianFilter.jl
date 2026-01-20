using OceananigansLagrangianFilter
using Oceananigans.Units
using Printf

# ### Model parameters
Nx = 10
Nz = 10
f = 1e-4                # Coriolis frequency [s⁻¹]
L_front = 10kilometers  # Initial front width [m]
aspect_ratio = 100      # L_front/H
Ro = 0.1                # Rossby number (defines M^2)

H = L_front/aspect_ratio  # Depth
M² = (Ro^2*f^2*L_front)/H # Horizontal buoyancy gradient
Δb = M²*L_front # Buoyancy difference across the front

filename_stem = "data/reference_sim";

# ### Define the grid
grid = RectilinearGrid(CPU(),size = (Nx, Nz), 
                       x = (-L_front/2, L_front/2),
                       z = (-H, 0),
                       topology = (Periodic, Flat, Bounded))

# ### Tracers
tracers = (:b,)

# ### Define the model
model =  NonhydrostaticModel(grid;
                coriolis = FPlane(f = f),
                buoyancy = BuoyancyTracer(),
                tracers = tracers,
                advection = WENO(),)


# ### Initialise the buoyancy and tracer 
bᵢ(x, z) = Δb*sin(2*pi/L_front * x)
set!(model, b= bᵢ) 

# ### Define the simulation
simulation = Simulation(model, Δt=20minutes, stop_time=1days)

# Set an adaptive timestep
conjure_time_step_wizard!(simulation, IterationInterval(20), cfl=0.2, max_Δt=20minutes)

# ### Set up the output 
u, v, w = model.velocities
b = model.tracers.b

# Output a jld2 file for Lagrangian filtering
simulation.output_writers[:jld2fields] = JLD2Writer(
    model, (; b, u, w), filename = filename_stem * ".jld2", schedule=TimeInterval(1hour), overwrite_existing=true)

# ### Run simulation
@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)


# ## Perform Lagrangian filtering
# Now we set up and run the offline Lagrangian filter on the output of the above simulation.


# ### Set up the filter configuration
filter_config = OfflineFilterConfig(original_data_filename="data/reference_sim.jld2", # Where the original simulation output is
                                    output_filename = "data/reference_offline_output.jld2", # Where to save the filtered output
                                    var_names_to_filter = ("b",), # Which variables to filter
                                    velocity_names = ("u","w"), # Velocities to use for Lagrangian filtering
                                    architecture = CPU(), # CPU() or GPU(), if GPU() make sure you have CUDA.jl installed and imported
                                    Δt = 20minutes, # Time step of filtering simulation
                                    T_out = 1hour, # How often to output filtered data
                                    N = 2, # Order of Butterworth filter
                                    freq_c = 1e-4/2, # Cut-off frequency of Butterworth filter
                                    compute_mean_velocities = true, # Whether to compute the mean velocities
                                    output_netcdf = false, # Whether to output filtered data to a netcdf file in addition to .jld2
                                    delete_intermediate_files = true, # Delete the individual output of the forward and backward passes
                                    compute_Eulerian_filter = true) # Whether to compute the Eulerian filter for comparison


# ### Run the offline Lagrangian filter
run_offline_Lagrangian_filter(filter_config)


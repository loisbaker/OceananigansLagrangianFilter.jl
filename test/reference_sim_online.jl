using OceananigansLagrangianFilter
using Oceananigans.Units

# Model parameters
Nx = 10
Nz = 10
f = 1e-4                # Coriolis frequency [s⁻¹]
L_front = 10kilometers  # Initial front width [m]
aspect_ratio = 100      # L_front/H
Ro = 0.1                # Rossby number (defines M^2)


H = L_front/aspect_ratio  # Depth
M² = (Ro^2*f^2*L_front)/H # Horizontal buoyancy gradient
Δb = M²*L_front # Buoyancy difference across the front


filename_stem = "data/reference_online_output";

# Grid
grid = RectilinearGrid(CPU(),size = (Nx, Nz), 
                       x = (-L_front/2, L_front/2),
                       z = (-H, 0),
                       topology = (Periodic, Flat, Bounded))

# Define model tracers
tracers = (:b,:T);

# Define model forcing
forcing = NamedTuple();

# Define filter configuration
filter_config = OnlineFilterConfig( grid = grid,
                                    output_filename = filename_stem * ".jld2",
                                    var_names_to_filter = ("b","T"), 
                                    velocity_names = ("u","w"),
                                    N = 1,
                                    freq_c = f/4)

# Create the filtered variables - these will be tracers in the model
filtered_vars = create_filtered_vars(filter_config)

# Add to the existing tracers
tracers = (filtered_vars..., tracers...)

# Create forcing for these filtered variables
filter_forcing = create_forcing(filtered_vars, filter_config);

# Add these to the existing forcing
forcing = merge(forcing, filter_forcing);

# Define closures
zero_filtered_var_closure = zero_closure_for_filtered_vars(filter_config)
horizontal_closure = HorizontalScalarDiffusivity(ν=1e-6, κ=merge((T=1e-6, b= 1e-6),zero_filtered_var_closure) )
vertical_closure = VerticalScalarDiffusivity(ν=1e-6 , κ=merge((T=1e-6, b= 1e-6),zero_filtered_var_closure) )
closure = (horizontal_closure, vertical_closure);

# Define the model
model =  NonhydrostaticModel(grid;
                coriolis = FPlane(f = f),
                buoyancy = BuoyancyTracer(),
                tracers = tracers,
                forcing = forcing,
                advection = WENO(),
                closure = closure)

# Initialise tracers
# Model buoyancy and tracers
bᵢ(x, z) = Δb*sin(2*pi/L_front * x)
Tᵢ(x, z) = exp(-(x/(L_front/50)).^2)
set!(model, b= bᵢ, T= Tᵢ) 

# Set appropriate initial conditions for the filtered variables based on the actual variables
initialise_filtered_vars_from_model(model, filter_config)

# ### Define the simulation
simulation = Simulation(model, Δt=20minutes, stop_time=1days)

# ### Set an adaptive timestep
conjure_time_step_wizard!(simulation, IterationInterval(20), cfl=0.2, max_Δt=20minutes)

# Create filtered outputs:
outputs = create_output_fields(model, filter_config);

# Add in original variables if needed:
outputs["b"] = model.tracers.b;
outputs["T"] = model.tracers.T;
outputs["u"] = model.velocities.u;
outputs["w"] = model.velocities.w;

# Output a .jld2 file:
simulation.output_writers[:jld2fields] = JLD2Writer(
    model, outputs, filename=filter_config.output_filename, schedule=TimeInterval(1hour), overwrite_existing=true)

# Run the simulation
@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)

# Option to regrid to mean position
regrid_to_mean_position!(filter_config)

# Option to calculate Eulerian filter too
compute_Eulerian_filter!(filter_config);

# Option to add a shifted time coordinate
compute_time_shift!(filter_config)

# Option to output a final NetCDF file
#jld2_to_netcdf(filename_stem * ".jld2", filename_stem * ".nc")



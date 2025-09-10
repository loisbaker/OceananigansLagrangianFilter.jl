
using Oceananigans
using Oceananigans.Models: ShallowWaterModel
using Printf
using NCDatasets
using CUDA

filename_stem = "SW_IO_with_tracer"

# Two-dimensional domain 
grid = RectilinearGrid(GPU(), size = (256, 256),
                       x = (0, 2*pi),
                       y = (0, 2*pi),
                       topology = (Periodic, Periodic, Flat))

# Building a `ShallowWaterModel`. We non-dimensionalise as in Kafiabad & Vanneste 2023.
Fr = 0.1 # Froude number
Ro = 1 # fRossby number

gravitational_acceleration = 1/Fr^2
coriolis = FPlane(f=1/Ro)

model = ShallowWaterModel(; grid, coriolis, gravitational_acceleration,
                            timestepper = :RungeKutta3,
                            tracers= (:T,),
                            momentum_advection = WENO())


# Set initial conditions 

# Velocity and height initial conditions - uniform velocity perturbation, initial height is 1 (unperturbed)
displacement = 2*pi/10
u_i = displacement/Ro   
h_i = 1
uh_i = u_i*h_i

# Initialise a tracer as a blob in the middle of the domain
width = 2*pi/15
T_i(x, y) = exp(-(((x - pi)^2 + (y - pi)^2)/width).^2)

set!(model, uh = uh_i, h= h_i, T = T_i )

## Build velocities 
uh, vh, h = model.solution

# Note that it is harder (and not implemented) to run a shallow water 
# simulation online, since the model evolves uh and vh rather than the velocities.

u = Field(uh / h)
v = Field(vh / h)
T = model.tracers.T

## Running a simulation
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

# Build the `output_writer` for the two-dimensional fields to be output.

# For Lagrangian filtering
simulation.output_writers[:fields_jld2] = JLD2Writer(model, (; u,v,T),
                                                        filename = filename_stem * ".jld2",
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)
# NetCDF can be useful too for visualisation
rm(filename_stem * ".nc",force=true)
simulation.output_writers[:fields_nc] = NetCDFWriter(model, (; u,v,T ),
                                                        filename = filename_stem * ".nc",
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)
#And finally run the simulation.
run!(simulation)

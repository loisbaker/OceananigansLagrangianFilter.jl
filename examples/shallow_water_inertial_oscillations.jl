
# Remember we've changed the source code in solution_and_tracer_tendencies.jl to give us modified shallow water

using Oceananigans
using Oceananigans.Models: ShallowWaterModel
using Printf
using NCDatasets
using CUDA

# Two-dimensional domain 

grid = RectilinearGrid(GPU(), size = (256, 256),
                       x = (0, 2*pi),
                       y = (0, 2*pi),
                       topology = (Periodic, Periodic, Flat))

# Building a `ShallowWaterModel`. We non-dimensionalise as in Kafiabaf & Vanneste 2023
Fr = 0.1
Ro = 1 

gravitational_acceleration = 1/Fr^2
coriolis = FPlane(f=1/Ro)

model = ShallowWaterModel(; grid, coriolis, gravitational_acceleration,
                            timestepper = :RungeKutta3,
                            tracers= (:T,),
                            momentum_advection = WENO())


f = coriolis.f
g = gravitational_acceleration

# Set initial conditions - uniform velocity perturbation, initial height is 1 (unperturbed)

max_displacement = 2*pi/10
u_i = max_displacement/Ro   
h_i = 1
uh_i = u_i*h_i

width = 2*pi/15
T_i(x, y) = exp(-((x - pi)/width).^2)

set!(model, uh = uh_i, h= h_i, T = T_i )

## Build velocities 
uh, vh, h = model.solution
u = Field(uh / h)
v = Field(vh / h)
T = model.tracers.T

## Running a simulation
simulation = Simulation(model, Δt = 1e-3, stop_time = 20)

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

fields_filename= joinpath(@__DIR__, "SW_IO_with_tracer")

# We'll output both jld2 and netcdf datasets
simulation.output_writers[:fields_jld2] = JLD2Writer(model, (; u,v,T),
                                                        filename = fields_filename,
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)

simulation.output_writers[:fields_nc] = NetCDFWriter(model, (; u,v,T ),
                                                        filename = fields_filename,
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)
#And finally run the simulation.
run!(simulation)


# Remember we've changed the source code in solution_and_tracer_tendencies.jl to give us modified shallow water

using Oceananigans
using Oceananigans.Models: ShallowWaterModel
using MAT
using Printf
using NCDatasets
using CUDA

# ## Two-dimensional domain 

# The shallow water model is two-dimensional and uses grids that are `Flat`
# in the vertical direction. 
arch = GPU()
grid = RectilinearGrid(arch, size = (256, 256),
                       x = (0, 2π),
                       y = (0, 2π),
                       topology = (Periodic, Periodic, Flat))

# Building a `ShallowWaterModel`

# We build a `ShallowWaterModel` with the `WENO` advection scheme,
# 3rd-order Runge-Kutta time-stepping, non-dimensional Coriolis, and
# gravitational acceleration

Fr = 0.3
Ro = 0.4 # These match geostrophic balance of the ICs, so don't change them without generating new ICs

gravitational_acceleration = 1/Fr^2
coriolis = FPlane(f=1/Ro)

model = ShallowWaterModel(; grid, coriolis, gravitational_acceleration,
                            timestepper = :RungeKutta3,
                            tracers= (:T,),
                            momentum_advection = WENO())


f = coriolis.f
g = gravitational_acceleration

# Get initial condition from mat file
data_file = joinpath(@__DIR__, "uvh_256_Fr_0_3_Ro_0_4_wave_0_5.mat")
file = matopen(data_file)
u_i = read(file, "u") 
v_i = read(file, "v")
h_i = read(file, "h")
close(file)


uh_i = u_i.*h_i 
vh_i = v_i.*h_i
width = 2*pi/15
T_i(x, y) = exp(-((x - pi)/width).^2)

set!(model, uh = uh_i, vh = vh_i, h= h_i, T = T_i )

## Build velocities and vorticity
uh, vh, h = model.solution
u = Field(uh / h)
v = Field(vh / h)
ω = Field(@at (Center, Center, Center) ∂x(v) - ∂y(u)) # We should output fields to filter at center locations
T = model.tracers.T

## Running a simulation
simulation = Simulation(model, Δt = 1e-3, stop_time = 40)

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

fields_filename= joinpath(@__DIR__, "SW_vort_with_tracer")

# We'll output both jld2 and netcdf datasets
simulation.output_writers[:fields_jld2] = JLD2Writer(model, (; ω,u,v,T),
                                                        filename = fields_filename,
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)

# simulation.output_writers[:fields_nc] = NetCDFWriter(model, (; ω,u,v,T ),
#                                                         filename = fields_filename,
#                                                         schedule = TimeInterval(0.1),
#                                                         overwrite_existing = true)
# And finally run the simulation.
run!(simulation)

#using CUDA: has_cuda
#using OrderedCollections: OrderedDict

using Oceananigans.Architectures: AbstractArchitecture
using Oceananigans.DistributedComputations: Distributed
using Oceananigans.Advection: Centered, adapt_advection_order
using Oceananigans.BoundaryConditions: regularize_field_boundary_conditions
using Oceananigans.Fields: Field, tracernames, VelocityFields, TracerFields, CenterField, location
using Oceananigans.Forcings: model_forcing
using Oceananigans.Grids: inflate_halo_size, with_halo, architecture
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Models: AbstractModel, NaNChecker, extract_boundary_conditions
using Oceananigans.TimeSteppers: Clock, TimeStepper, update_state!
using Oceananigans.TurbulenceClosures: validate_closure, with_tracers, build_diffusivity_fields, time_discretization, implicit_diffusion_solver
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: FlavorOfCATKE
using Oceananigans.Utils: tupleit
using Oceananigans.Grids: topology
import Oceananigans.Architectures: architecture
import Oceananigans.Models: total_velocities, default_nan_checker, timestepper



mutable struct LagrangianFilter{TS, E, A<:AbstractArchitecture, G, T, U, C, F,
                                   V, K, B, AF} <: AbstractModel{TS, A}

         architecture :: A        # Computer `Architecture` on which `Model` is run
                 grid :: G        # Grid of physical points on which `Model` is solved
                clock :: Clock{T} # Tracks iteration number and simulation time of `Model`
            advection :: V        # Advection scheme for velocities _and_ tracers
              forcing :: F        # Container for forcing functions defined by the user
              closure :: E        # Diffusive 'turbulence closure' for all model fields
            buoyancy  :: B        # Defunct buoyancy
           velocities :: U        # Container for velocity fields `u`, `v`, and `w`
              tracers :: C        # Container for tracer fields
   diffusivity_fields :: K        # Container for turbulent diffusivities
          timestepper :: TS       # Object containing timestepper fields and parameters
     auxiliary_fields :: AF       # User-specified auxiliary fields for forcing functions and boundary conditions
end

"""
    LagrangianFilter(;           grid,
                                    clock = Clock{eltype(grid)}(time = 0),
                                advection = Centered(),
                      forcing::NamedTuple = NamedTuple(),
                                  closure = nothing,
          boundary_conditions::NamedTuple = NamedTuple(),
                                  tracers = (),
                              timestepper = :RungeKutta3,
                               velocities = nothing,
                       diffusivity_fields = nothing,
                         auxiliary_fields = NamedTuple())

By default, all Bounded directions are rigid and impenetrable.

Keyword arguments
=================

  - `grid`: (required) The resolution and discrete geometry on which the `model` is solved. The
            architecture (CPU/GPU) that the model is solved on is inferred from the architecture
            of the `grid`. Note that the grid needs to be regularly spaced in the horizontal
            dimensions, ``x`` and ``y``.
  - `advection`: The scheme that advects velocities and tracers. See `Oceananigans.Advection`.
  - `forcing`: `NamedTuple` of user-defined forcing functions that contribute to solution tendencies.
  - `closure`: The turbulence closure for `model`. See `Oceananigans.TurbulenceClosures`.
  - `boundary_conditions`: `NamedTuple` containing field boundary conditions.
  - `tracers`: A tuple of symbols defining the names of the modeled tracers, or a `NamedTuple` of
               preallocated `CenterField`s.
  - `timestepper`: A symbol that specifies the time-stepping method. Either `:QuasiAdamsBashforth2` or
                   `:RungeKutta3` (default).
  - `velocities`: The model velocities. Default: `nothing`.
  - `diffusivity_fields`: Diffusivity fields. Default: `nothing`.
  - `auxiliary_fields`: `NamedTuple` of auxiliary fields. Default: `nothing`         
"""


function LagrangianFilter(; grid,
                            clock = Clock{eltype(grid)}(time = 0),
                            advection = Centered(),
                            forcing::NamedTuple = NamedTuple(),
                            closure = nothing,
                            boundary_conditions::NamedTuple = NamedTuple(),
                            tracers = (),
                            timestepper = :RungeKutta3,
                            velocities = nothing,
                            diffusivity_fields = nothing,
                            auxiliary_fields = NamedTuple())

    arch = architecture(grid)

    tracers = tupleit(tracers) # supports tracers=:c keyword argument (for example)

   
    # We don't support CAKTE for LagrangianFilter yet.
    closure = validate_closure(closure)
    first_closure = closure isa Tuple ? first(closure) : closure
    first_closure isa FlavorOfCATKE &&
        error("CATKEVerticalDiffusivity is not supported for LagrangianFilter --- yet!")

    # Adjust advection scheme to be valid on a particular grid size. i.e. if the grid size
    # is smaller than the advection order, reduce the order of the advection in that particular
    # direction
    advection = adapt_advection_order(advection, grid)

    # Adjust halos when the advection scheme or turbulence closure requires it.
    # Note that halos are isotropic by default; however we respect user-input here
    # by adjusting each (x, y, z) halo individually.
    grid = inflate_grid_halo_size(grid, advection, closure)

    # Collect boundary conditions for all model prognostic fields and, if specified, some model
    # auxiliary fields. Boundary conditions are "regularized" based on the _name_ of the field:
    # boundary conditions on u, v, w are regularized assuming they represent momentum at appropriate
    # staggered locations. All other fields are regularized assuming they are tracers.
    # Note that we do not regularize boundary conditions contained in *tupled* diffusivity fields right now.

    # First, we extract boundary conditions that are embedded within any _user-specified_ field tuples:
    embedded_boundary_conditions = merge(extract_boundary_conditions(velocities),
                                         extract_boundary_conditions(tracers),
                                         extract_boundary_conditions(diffusivity_fields))

    # Next, we form a list of default boundary conditions:
    prognostic_field_names = (:u, :v, :w, tracernames(tracers)..., keys(auxiliary_fields)...)
    default_boundary_conditions = NamedTuple{prognostic_field_names}(FieldBoundaryConditions()
                                                                     for name in prognostic_field_names)

    # Finally, we merge specified, embedded, and default boundary conditions. Specified boundary conditions
    # have precedence, followed by embedded, followed by default.
    boundary_conditions = merge(default_boundary_conditions, embedded_boundary_conditions, boundary_conditions)
    boundary_conditions = regularize_field_boundary_conditions(boundary_conditions, grid, prognostic_field_names)

    # Ensure `closure` describes all tracers
    closure = with_tracers(tracernames(tracers), closure)

    # If grid has been inflated (e.g. to accommodate different advection scheme to original data), the auxiliary_fields
    # are not automatically updated to reflect this, so we do this here
    auxiliary_fields_update = NamedTuple()
    for (name, field) in pairs(auxiliary_fields)
        loc = location(field)
        new_field = Field((loc[1](), loc[2](), loc[3]()), grid)
        set!(new_field, field)
        auxiliary_fields_update = merge(auxiliary_fields_update, (;(name => new_field)))
    end        
    auxiliary_fields = auxiliary_fields_update     

    # Either check grid-correctness, or construct tuples of fields
    velocities         = VelocityFields(velocities, grid, boundary_conditions)
    tracers            = TracerFields(tracers,      grid, boundary_conditions)
    diffusivity_fields = build_diffusivity_fields(diffusivity_fields, grid, clock, tracernames(tracers), boundary_conditions, closure)
    buoyancy = nothing                                                                    
    model_fields = merge(velocities, tracers, auxiliary_fields)
    prognostic_fields = merge(velocities, tracers)

    # Instantiate timestepper if not already instantiated
    implicit_solver = implicit_diffusion_solver(time_discretization(closure), grid)
    timestepper = TimeStepper(timestepper, grid, prognostic_fields; implicit_solver=implicit_solver)

    # Regularize forcing for model tracer and velocity fields.
    forcing = model_forcing(model_fields; forcing...)

    model = LagrangianFilter(arch, grid, clock, advection,
                              forcing, closure, buoyancy, velocities, tracers,
                              diffusivity_fields, timestepper, auxiliary_fields)

    update_state!(model; compute_tendencies = false)
    
    return model
end

architecture(model::LagrangianFilter) = model.architecture

function inflate_grid_halo_size(grid, tendency_terms...)
    user_halo = grid.Hx, grid.Hy, grid.Hz
    required_halo = Hx, Hy, Hz = inflate_halo_size(user_halo..., grid, tendency_terms...)

    if any(user_halo .< required_halo) # Replace grid
        @warn "Inflating model grid halo size to ($Hx, $Hy, $Hz) and recreating grid. " *
              "Note that an ImmersedBoundaryGrid requires an extra halo point in all non-flat directions compared to a non-immersed boundary grid."
              "The model grid will be different from the input grid. To avoid this warning, " *
              "pass halo=($Hx, $Hy, $Hz) when constructing the grid."

        grid = with_halo((Hx, Hy, Hz), grid)
    end

    return grid
end

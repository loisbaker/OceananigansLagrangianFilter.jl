# Regridding the filtered fields to the mean position.
#
# Each filtered field f* is stored at the grid points x, which label the particles there, but belongs at the particles'
# mean positions X(x) = x + xi(x). For each grid point y:
#
#   1. find the label x* with X(x*) = y, by Newton's method, with the maps xi interpolated linearly between cell centres;
#   2. if x* is in the fluid, the regridded field at y is f*(x*), interpolated in the same way;
#   3. otherwise no particle in the fluid has its mean position at y, and the regridded field there is NaN.
#
# Labels are fractional cell-centre indices (label i is cell centre i), so stretched grids only enter through the
# coordinates of the cell centres. Every array has one ghost cell beyond each end of every non-Flat dimension: wrapped
# in a periodic dimension, and extrapolated linearly beyond a bounded boundary, so nothing is assumed about the flow at
# the boundary. A label is in the fluid if it is inside the domain, in a fluid cell, and above the bottom (for an
# immersed boundary with a bottom height). The same kernel runs on CPU and GPU.

const FLAT_DIMENSION, PERIODIC_DIMENSION, BOUNDED_DIMENSION = 0, 1, 2

# Status of each regridded point. Any status but REGRID_OK gives NaN.
const REGRID_OK, REGRID_SOLID_POINT, REGRID_NO_MEAN_POSITION, REGRID_NOT_CONVERGED = 0, 1, 2, 3

# The architecture to regrid on: the offline filter's, or the online model grid's
regrid_architecture(config::AbstractOfflineConfig) = config.architecture
regrid_architecture(config::AbstractOnlineConfig) = Oceananigans.Architectures.architecture(config.grid)

#####
##### Kernel helpers. Arrays are padded with one ghost cell beyond each end of every non-Flat dimension, so cell
##### centre i is at padded index i + 1 (or 1 in a Flat dimension).
#####

# The padded index of the lower neighbour of the fractional label f in one dimension, the weight of the upper
# neighbour, and the offset to the upper neighbour (zero in a Flat dimension). Periodic labels are wrapped, bounded
# labels are kept within the ghost cells.
@inline function label_neighbours(f, N, kind)
    f = ifelse(kind == PERIODIC_DIMENSION, mod(f - 1/2, N) + 1/2, clamp(f, zero(f), N + one(f)))
    i = min(floor(Int, f), N)
    flat = kind == FLAT_DIMENSION
    return (ifelse(flat, 1, i + 1), ifelse(flat, zero(f), f - i), ifelse(flat, 0, 1))
end

# Trilinear interpolation of the padded array a (with trailing indices m) at the neighbours c, and its derivatives
# with respect to the three labels
@inline function trilinear_with_gradient(a, c, m...)
    (i, wi, si), (j, wj, sj), (k, wk, sk) = c
    @inbounds begin
        a000 = a[i, j, k, m...];       a100 = a[i+si, j, k, m...]
        a010 = a[i, j+sj, k, m...];    a110 = a[i+si, j+sj, k, m...]
        a001 = a[i, j, k+sk, m...];    a101 = a[i+si, j, k+sk, m...]
        a011 = a[i, j+sj, k+sk, m...]; a111 = a[i+si, j+sj, k+sk, m...]
    end
    a00 = a000 + wi * (a100 - a000); a10 = a010 + wi * (a110 - a010)
    a01 = a001 + wi * (a101 - a001); a11 = a011 + wi * (a111 - a011)
    a0 = a00 + wj * (a10 - a00);     a1 = a01 + wj * (a11 - a01)
    dx0 = (a100 - a000) + wj * ((a110 - a010) - (a100 - a000))
    dx1 = (a101 - a001) + wj * ((a111 - a011) - (a101 - a001))
    dx = dx0 + wk * (dx1 - dx0)
    dy = (a10 - a00) + wk * ((a11 - a01) - (a10 - a00))
    dz = a1 - a0
    return a0 + wk * (a1 - a0), (dx, dy, dz)
end

@inline trilinear(a, c, m...) = trilinear_with_gradient(a, c, m...)[1]

# One component of the mean position X(f) - y (wrapped into the nearest period), and its derivatives with respect to
# the three labels (row d of the Jacobian). Rows of Flat dimensions are the identity, so they never move.
@inline function residual_row(xi, nodes, c, y, kind, period, d)
    x, (gx, gy, gz) = trilinear_with_gradient(xi, c)
    i, w, s = c[d]
    @inbounds spacing = nodes[i+s] - nodes[i]
    @inbounds r = nodes[i] + w * spacing + x - y
    r = ifelse(kind == PERIODIC_DIMENSION, r - period * round(r / period), r)
    flat = kind == FLAT_DIMENSION
    o, z = one(r), zero(r)
    J = (ifelse(flat, ifelse(d == 1, o, z), gx + ifelse(d == 1, spacing, z)),
         ifelse(flat, ifelse(d == 2, o, z), gy + ifelse(d == 2, spacing, z)),
         ifelse(flat, ifelse(d == 3, o, z), gz + ifelse(d == 3, spacing, z)))
    return ifelse(flat, z, r), J
end

@inline function mean_position_residual(f, y, xi, nodes, kinds, Ns, periods)
    c = (label_neighbours(f[1], Ns[1], kinds[1]), label_neighbours(f[2], Ns[2], kinds[2]),
         label_neighbours(f[3], Ns[3], kinds[3]))
    r1, J1 = residual_row(xi[1], nodes[1], c, y[1], kinds[1], periods[1], 1)
    r2, J2 = residual_row(xi[2], nodes[2], c, y[2], kinds[2], periods[2], 2)
    r3, J3 = residual_row(xi[3], nodes[3], c, y[3], kinds[3], periods[3], 3)
    return (r1, r2, r3), (J1, J2, J3)
end

@inline det3(a, b, c) = a[1] * (b[2] * c[3] - b[3] * c[2]) - a[2] * (b[1] * c[3] - b[3] * c[1]) + a[3] * (b[1] * c[2] - b[2] * c[1])

# The Newton step -J \ r by Cramer's rule, or a fixed-point step -r / h if J is (nearly) singular
@inline function newton_step(J, r, h)
    c1 = (J[1][1], J[2][1], J[3][1]); c2 = (J[1][2], J[2][2], J[3][2]); c3 = (J[1][3], J[2][3], J[3][3])
    D = det3(c1, c2, c3)
    singular = abs(D) < 1e-6 * abs(h[1] * h[2] * h[3])
    return (ifelse(singular, -r[1] / h[1], -det3(r, c2, c3) / D),
            ifelse(singular, -r[2] / h[2], -det3(c1, r, c3) / D),
            ifelse(singular, -r[3] / h[3], -det3(c1, c2, r) / D))
end

# The squared residual in cells
@inline merit(r, h) = (r[1] / h[1])^2 + (r[2] / h[2])^2 + (r[3] / h[3])^2

@inline padding(kind) = ifelse(kind == FLAT_DIMENSION, 0, 1)

# The grid spacing at cell centre i (padded index), or 1 in a Flat dimension
@inline centre_spacing(nodes, i, kind) = ifelse(kind == FLAT_DIMENSION, one(eltype(nodes)),
                                                @inbounds((nodes[i+padding(kind)] - nodes[i-padding(kind)]) / 2))

# Wrap a periodic label, and keep a bounded label within the ghost cells, as in label_neighbours
@inline tidy_label(f, N, kind) = ifelse(kind == PERIODIC_DIMENSION, mod(f - 1/2, N) + 1/2,
                                        ifelse(kind == BOUNDED_DIMENSION, clamp(f, zero(f), N + one(f)), f))

@inline inside_domain(f, N, kind, tol) = (kind != BOUNDED_DIMENSION) | ((f >= 1/2 - tol) & (f <= N + 1/2 + tol))

# The padded index of the cell containing the label f
@inline nearest_centre(f, N, kind) = clamp(round(Int, f), 1, N) + padding(kind)

# Whether the label f (with neighbours c) is in the fluid: inside the domain, in a fluid cell, and above the bottom
@inline function label_in_fluid(f, c, nodes, wet, bottom, kinds, Ns, tol)
    inside = inside_domain(f[1], Ns[1], kinds[1], tol) & inside_domain(f[2], Ns[2], kinds[2], tol) &
             inside_domain(f[3], Ns[3], kinds[3], tol)
    in_fluid_cell = @inbounds wet[nearest_centre(f[1], Ns[1], kinds[1]), nearest_centre(f[2], Ns[2], kinds[2]),
                                  nearest_centre(f[3], Ns[3], kinds[3])]
    k, w, s = c[3]
    @inbounds z = nodes[3][k] + w * (nodes[3][k+s] - nodes[3][k])
    above_bottom = z >= trilinear(bottom, (c[1], c[2], (1, zero(z), 0)))
    return inside & in_fluid_cell & above_bottom
end

# Newton's method with backtracking (Armijo) for the label x* with X(x*) = y, then the fields interpolated at x*.
# The residual is measured in cells, using the grid spacing h at y.
@kernel function regrid_to_mean_kernel!(out, status, fields, xi, nodes, wet, bottom, kinds, Ns, periods, maxiter, tol)
    i, j, k = @index(Global, NTuple)
    FT = eltype(out)
    ip, jp, kp = i + padding(kinds[1]), j + padding(kinds[2]), k + padding(kinds[3])
    @inbounds y = (nodes[1][ip], nodes[2][jp], nodes[3][kp])
    h = (centre_spacing(nodes[1], ip, kinds[1]), centre_spacing(nodes[2], jp, kinds[2]), centre_spacing(nodes[3], kp, kinds[3]))

    # 1. Newton's method for the label x* with mean position y, if y is in the fluid
    f = (FT(i), FT(j), FT(k))
    fluid_point = @inbounds wet[ip, jp, kp] & (y[3] >= bottom[ip, jp, 1])
    converged = true
    if fluid_point
        r, J = mean_position_residual(f, y, xi, nodes, kinds, Ns, periods)
        m = merit(r, h)
        n = 0
        while (m > tol^2) & (n < maxiter)
            n += 1
            df = newton_step(J, r, h)
            t = one(FT)
            new_f, new_r, new_J, new_m = f, r, J, m
            for _ in 1:20 # halve the step until the residual falls enough
                new_f = (f[1] + t * df[1], f[2] + t * df[2], f[3] + t * df[3])
                new_r, new_J = mean_position_residual(new_f, y, xi, nodes, kinds, Ns, periods)
                new_m = merit(new_r, h)
                new_m <= (1 - FT(1e-4) * t) * m && break
                t /= 2
            end
            f, r, J, m = new_f, new_r, new_J, new_m
        end
        converged = m <= tol^2
    end
    f = (tidy_label(f[1], Ns[1], kinds[1]), tidy_label(f[2], Ns[2], kinds[2]), tidy_label(f[3], Ns[3], kinds[3]))
    c = (label_neighbours(f[1], Ns[1], kinds[1]), label_neighbours(f[2], Ns[2], kinds[2]),
         label_neighbours(f[3], Ns[3], kinds[3]))

    # 2. Is x* in the fluid? (Test this before convergence: a label outside the fluid may not be reachable, since
    #    labels stop at the ghost cells.)
    s = ifelse(!fluid_point, REGRID_SOLID_POINT,
        ifelse(!label_in_fluid(f, c, nodes, wet, bottom, kinds, Ns, tol), REGRID_NO_MEAN_POSITION,
        ifelse(!converged, REGRID_NOT_CONVERGED, REGRID_OK)))
    @inbounds status[i, j, k] = s

    # 3. The fields at x*, or NaN
    for v in 1:size(fields, 4)
        @inbounds out[i, j, k, v] = ifelse(s == REGRID_OK, trilinear(fields, c, v), FT(NaN))
    end
end

# The mean position of the label at padded indices (i, j, k)
@inline mean_position(xi, nodes, i, j, k) =
    @inbounds (nodes[1][i] + xi[1][i, j, k], nodes[2][j] + xi[2][i, j, k], nodes[3][k] + xi[3][i, j, k])

# The edge of a deformed cell along dimension d from the corner (i, j, k), oriented along increasing labels, where
# (i0, j0, k0) is the cell's lower corner and s the offsets to its upper corner. In a Flat dimension it is a unit vector.
@inline function cell_edge(xi, nodes, i, j, k, i0, j0, k0, s, d)
    FT = eltype(nodes[1])
    lower = mean_position(xi, nodes, ifelse(d == 1, i0, i), ifelse(d == 2, j0, j), ifelse(d == 3, k0, k))
    upper = mean_position(xi, nodes, ifelse(d == 1, i0 + s[1], i), ifelse(d == 2, j0 + s[2], j), ifelse(d == 3, k0 + s[3], k))
    flat = s[d] == 0
    o, z = one(FT), zero(FT)
    return (ifelse(flat, ifelse(d == 1, o, z), upper[1] - lower[1]),
            ifelse(flat, ifelse(d == 2, o, z), upper[2] - lower[2]),
            ifelse(flat, ifelse(d == 3, o, z), upper[3] - lower[3]))
end

# The smallest Jacobian determinant of the mean-position map at the corners of the deformed cell between label
# centres i:i+1, j:j+1, k:k+1 (Inf for cells that cross a bounded boundary or have a corner label outside the fluid,
# as in label_in_fluid). The map folds where this is <= 0.
@kernel function folded_cells_kernel!(min_det, xi, nodes, wet, bottom, kinds, Ns)
    i, j, k = @index(Global, NTuple)
    FT = eltype(min_det)
    s = (padding(kinds[1]), padding(kinds[2]), padding(kinds[3])) # the offsets to the upper corner
    i0, j0, k0 = i + s[1], j + s[2], k + s[3]
    valid = ((kinds[1] != BOUNDED_DIMENSION) | (i < Ns[1])) & ((kinds[2] != BOUNDED_DIMENSION) | (j < Ns[2])) &
            ((kinds[3] != BOUNDED_DIMENSION) | (k < Ns[3]))
    all_fluid = true
    D = FT(Inf)
    for a in (0, s[1]), b in (0, s[2]), c in (0, s[3])
        ic, jc, kc = i0 + a, j0 + b, k0 + c
        all_fluid &= @inbounds wet[ic, jc, kc] & (nodes[3][kc] >= bottom[ic, jc, 1])
        D = min(D, det3(cell_edge(xi, nodes, ic, jc, kc, i0, j0, k0, s, 1),
                        cell_edge(xi, nodes, ic, jc, kc, i0, j0, k0, s, 2),
                        cell_edge(xi, nodes, ic, jc, kc, i0, j0, k0, s, 3)))
    end
    @inbounds min_det[i, j, k] = ifelse(valid & all_fluid, D, FT(Inf))
end

#####
##### Setting up on the host: grid geometry, ghost cells and the maps
#####

dimension_kind(::Type{Flat}) = FLAT_DIMENSION
dimension_kind(::Type{Periodic}) = PERIODIC_DIMENSION
dimension_kind(::Type{Bounded}) = BOUNDED_DIMENSION

# The coordinates along dimension d at location loc: longitude and latitude in degrees on a LatitudeLongitudeGrid
# (whose xnodes and ynodes are distances in metres)
regrid_coordinates(grid::LatitudeLongitudeGrid, d, loc) = collect((λnodes, φnodes, znodes)[d](grid, loc))
regrid_coordinates(grid, d, loc) = collect((xnodes, ynodes, znodes)[d](grid, loc))

# Cell-centre coordinates, with a ghost centre beyond each end: one period away for a periodic dimension, mirrored in
# the boundary for a bounded one. Longitude and latitude are in degrees.
function padded_nodes(grid, d, kind)
    grid = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    kind == FLAT_DIMENSION && return [zero(eltype(grid))]
    centres = regrid_coordinates(grid, d, Center())
    faces = regrid_coordinates(grid, d, Face()) # only N faces in a periodic dimension
    N = length(centres)
    if kind == PERIODIC_DIMENSION
        L = (grid.Lx, grid.Ly, grid.Lz)[d]
        return [centres[N] - L; centres; centres[1] + L]
    else
        return [2faces[1] - centres[1]; centres; 2faces[end] - centres[N]]
    end
end

# Fill the ghost cells of the padded array a along dimension d: wrapped for a periodic dimension, and extrapolated
# linearly beyond each bounded boundary
function fill_ghosts!(a, d, kind, nodes)
    kind == FLAT_DIMENSION && return a
    N = size(a, d) - 2
    ghost(i) = selectdim(a, d, i)
    if kind == PERIODIC_DIMENSION
        ghost(1) .= ghost(N + 1)
        ghost(N + 2) .= ghost(2)
    elseif N == 1
        ghost(1) .= ghost(2)
        ghost(3) .= ghost(2)
    else
        for (g, i1, i2) in ((1, 2, 3), (N + 2, N + 1, N))
            ghost(g) .= ghost(i1) .+ (ghost(i1) .- ghost(i2)) .* ((nodes[i1] - nodes[g]) / (nodes[i2] - nodes[i1]))
        end
    end
    return a
end

# Give solid cells next to the fluid the mean of their fluid (or already filled) neighbours, two cells deep, so that
# interpolation near an immersed boundary only uses fluid values. Deeper solid cells are never used.
function fill_solid_cells!(a, wet)
    filled = copy(wet)
    for _ in 1:2
        new_filled = copy(filled)
        for I in CartesianIndices(wet)
            filled[I] && continue
            total, count = zero(eltype(a)), 0
            for d in 1:3, step in (-1, 1)
                J = I + CartesianIndex(ntuple(e -> e == d ? step : 0, 3))
                checkbounds(Bool, wet, J) && filled[J] || continue
                total += a[J]; count += 1
            end
            count > 0 && (a[I] = total / count; new_filled[I] = true)
        end
        filled = new_filled
    end
    return a
end

# The interior of an array with halos, padded with ghost cells
function padded_interior(data, grid, setup)
    interior = _remove_halos(data, grid)
    setup.immersed && fill_solid_cells!(interior, setup.wet_interior)
    a = zeros(eltype(interior), size(interior) .+ 2 .* setup.pads)
    a[ntuple(d -> setup.pads[d] .+ (1:size(interior, d)), 3)...] .= interior
    for d in 1:3
        fill_ghosts!(a, d, setup.kinds[d], setup.nodes[d])
    end
    return a
end

function regrid_setup(grid)
    kinds = map(dimension_kind, topology(grid))
    Ns = size(grid)
    pads = map(kind -> kind == FLAT_DIMENSION ? 0 : 1, kinds)
    nodes = ntuple(d -> padded_nodes(grid, d, kinds[d]), 3)
    periods = ntuple(d -> kinds[d] == PERIODIC_DIMENSION ? nodes[d][end] - nodes[d][2] : one(eltype(grid)), 3) # 1 is unused
    immersed = grid isa ImmersedBoundaryGrid
    wet_interior = immersed ? [!immersed_cell(i, j, k, grid) for i in 1:Ns[1], j in 1:Ns[2], k in 1:Ns[3]] : trues(Ns)
    wet = falses(Ns .+ 2 .* pads)
    wet[ntuple(d -> pads[d] .+ (1:Ns[d]), 3)...] .= wet_interior
    for d in 1:3 # ghost cells are wet if the cell they copy or extrapolate from is
        kinds[d] == FLAT_DIMENSION && continue
        N = Ns[d]
        lo, hi = kinds[d] == PERIODIC_DIMENSION ? (N + 1, 2) : (2, N + 1)
        selectdim(wet, d, 1) .= selectdim(wet, d, lo); selectdim(wet, d, N + 2) .= selectdim(wet, d, hi)
    end
    # The bottom height in each column, padded in x and y (the lowest finite number without an immersed bottom, which
    # unlike -Inf can be interpolated)
    FT = eltype(grid)
    bottom = fill(-floatmax(FT), Ns[1] + 2pads[1], Ns[2] + 2pads[2], 1)
    if immersed && hasproperty(grid.immersed_boundary, :bottom_height)
        bottom[pads[1] .+ (1:Ns[1]), pads[2] .+ (1:Ns[2]), 1] .= [grid.immersed_boundary.bottom_height[i, j, 1] for i in 1:Ns[1], j in 1:Ns[2]]
        fill_ghosts!(bottom, 1, kinds[1], nodes[1])
        fill_ghosts!(bottom, 2, kinds[2], nodes[2])
    end
    return (; kinds, Ns, pads, nodes, periods, immersed, wet_interior, wet, bottom)
end

# The maps (xi_u, xi_v, xi_w), given with halos (or `nothing` for no map), padded, in the units of each coordinate
# (degrees for longitude and latitude). A dimension without a map (or Flat) has zero displacement.
function padded_maps(maps, grid, setup)
    FT = eltype(grid)
    underlying = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    return ntuple(3) do d
        (setup.kinds[d] == FLAT_DIMENSION || isnothing(maps[d])) && return zeros(FT, size(setup.wet))
        data = FT.(maps[d])
        if underlying isa LatitudeLongitudeGrid && d < 3 # metres to degrees, at each cell's latitude
            R = underlying.radius
            data = d == 1 ? data ./ (R .* cosd.(reshape(collect(φnodes(underlying, Center(), with_halos = true)), 1, :, 1))) .* (180 / pi) :
                            data ./ R .* (180 / pi)
        end
        padded_interior(data, grid, setup)
    end
end

# Regrid the fields at one output time. `maps` are the three maps (xi_u, xi_v, xi_w) and `fields` the filtered fields,
# all with halos as saved in the output file. Returns the regridded fields (interior only, the last index numbering the
# fields), the status of each point, and the number of folded cells.
function regrid_fields(grid, maps, fields; architecture = CPU(), maxiter = 50, setup = regrid_setup(grid))
    FT = eltype(grid)
    Ns = setup.Ns
    nodes = map(n -> on_architecture(architecture, FT.(n)), setup.nodes)
    wet = on_architecture(architecture, Array{Bool}(setup.wet))
    bottom = on_architecture(architecture, setup.bottom)
    xi = map(a -> on_architecture(architecture, a), padded_maps(maps, grid, setup))
    padded_fields = on_architecture(architecture, cat((padded_interior(FT.(f), grid, setup) for f in fields)...; dims = 4))
    tol = sqrt(eps(FT)) # in cells

    backend = Oceananigans.Architectures.device(architecture)
    regridded = on_architecture(architecture, zeros(FT, Ns..., length(fields)))
    status = on_architecture(architecture, zeros(Int8, Ns))
    regrid_to_mean_kernel!(backend)(regridded, status, padded_fields, xi, nodes, wet, bottom, setup.kinds, Ns, setup.periods,
                                    maxiter, tol; ndrange = Ns)
    min_det = on_architecture(architecture, zeros(FT, Ns))
    folded_cells_kernel!(backend)(min_det, xi, nodes, wet, bottom, setup.kinds, Ns; ndrange = Ns)
    return Array(regridded), Array(status), count(<=(0), min_det)
end

"""
    regrid_to_mean_position!(config::AbstractConfig; extra_vars_to_regrid = (), architecture = regrid_architecture(config),
                             maxiter = 50)

Regrid the Lagrangian filtered fields to the mean position, and save each as a new variable with the suffix
`_at_mean` in the output file.

Each filtered field ``f^*`` is stored at the regular grid points, which label the particles, but belongs at their mean
position ``\\Xi(\\boldsymbol{x}) = \\boldsymbol{x} + \\boldsymbol{\\xi}(\\boldsymbol{x})``. For each grid point
``\\boldsymbol{y}``, the label ``\\boldsymbol{x}^*`` with ``\\Xi(\\boldsymbol{x}^*) = \\boldsymbol{y}`` is found by
Newton's method, with the maps ``\\boldsymbol{\\xi}`` interpolated linearly between cell centres, and the regridded
field is ``f^*(\\boldsymbol{x}^*)``. This runs on the CPU or GPU.

The regridded field is `NaN` at grid points that are not the mean position of any particle in the fluid: where the
label ``\\boldsymbol{x}^*`` lies outside the domain, or in an immersed solid (in a solid cell, or below the bottom
height). This happens, for example, near a boundary where the particles' mean positions are pulled away from it. It is
also `NaN` at solid grid points and where Newton's method does not converge. Beyond bounded boundaries, the maps and
fields are extrapolated linearly, so nothing is assumed about the flow at the boundary. A warning is given if the map
folds (so that some grid points are the mean position of more than one label).

Arguments
=========
- `config`: The `OfflineFilterConfig` or `OnlineFilterConfig` used for the filtering.
- `extra_vars_to_regrid`: Additional filtered variables in the output file to regrid.
- `architecture`: The architecture to regrid on. Default: the offline filter's `architecture`, or the architecture of
  the online model's grid.
- `maxiter`: The maximum number of Newton iterations. Default: `50`.
"""
function regrid_to_mean_position!(config::AbstractConfig; extra_vars_to_regrid::Tuple{Vararg{String}} = (),
                                  architecture = regrid_architecture(config), maxiter::Int = 50)
    config.lagrangian || error("Regridding to mean position is not meaningful for Eulerian filtering")
    var_names_to_filter = config.compute_mean_velocities ? (config.var_names_to_filter..., config.velocity_names...) :
                                                           config.var_names_to_filter
    var_names_to_regrid = Tuple(unique((Tuple(var * config.label * "_Lagrangian_filtered" for var in var_names_to_filter)...,
                                        extra_vars_to_regrid...)))
    output_filename = config.output_filename
    counts = zeros(Int, 3)
    folded_cells = 0
    folded_times = 0

    jldopen(output_filename, "r+") do file
        grid = file["serialized/grid"]
        underlying = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
        underlying isa Union{RectilinearGrid, LatitudeLongitudeGrid} ||
            error("Regridding to mean position is only implemented for RectilinearGrid and LatitudeLongitudeGrid")
        for var in var_names_to_regrid
            file["timeseries/$var/serialized/location"] == (Center, Center, Center) ||
                error("Regridding to mean position needs $var at cell centres")
        end
        FT = eltype(grid)
        setup = regrid_setup(grid)

        # Copy the serialized properties of each variable for its regridded version
        for var in var_names_to_regrid
            new_path = "timeseries/$(var)_at_mean/serialized"
            haskey(file, new_path) && Base.delete!(file, new_path)
            g = Group(file, new_path)
            for property in keys(file["timeseries/$var/serialized"])
                g[property] = file["timeseries/$var/serialized/$property"]
            end
        end

        Ns = setup.Ns
        halo = (grid.Hx, grid.Hy, grid.Hz)
        for iter in parse.(Int, keys(file["timeseries/t"]))
            maps = ntuple(3) do d
                vel = ("u", "v", "w")[d]
                vel in config.velocity_names ? file["timeseries/xi_$vel$(config.label)/$iter"] : nothing
            end
            fields = Tuple(file["timeseries/$var/$iter"] for var in var_names_to_regrid)
            out, status, n_folded = regrid_fields(grid, maps, fields; architecture, maxiter, setup)
            for s in 1:3
                counts[s] += count(==(s), status)
            end
            folded_cells += n_folded
            folded_times += n_folded > 0

            for (v, var) in enumerate(var_names_to_regrid)
                data = zeros(FT, Ns .+ 2 .* halo)
                data[ntuple(d -> halo[d] .+ (1:Ns[d]), 3)...] .= out[:, :, :, v]
                new_var_loc = "timeseries/$(var)_at_mean/$iter"
                haskey(file, new_var_loc) && Base.delete!(file, new_var_loc)
                file[new_var_loc] = data
            end
        end
    end

    counts[REGRID_NOT_CONVERGED] > 0 &&
        @warn "Regridding to mean position: Newton's method did not converge at $(counts[REGRID_NOT_CONVERGED]) points (set to NaN)"
    counts[REGRID_NO_MEAN_POSITION] > 0 &&
        @info "Regridding to mean position: $(counts[REGRID_NO_MEAN_POSITION]) points are not the mean position of any particle in the fluid (set to NaN)"
    folded_cells > 0 &&
        @warn "Regridding to mean position: the map folds in $folded_cells cells over $folded_times output times, so some points are the mean position of more than one label there"
    @info "Wrote regridded data to new variables with _at_mean suffix in file $output_filename"
    return nothing
end

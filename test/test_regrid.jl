using OceananigansLagrangianFilter.Utils: regrid_fields, REGRID_OK, REGRID_SOLID_POINT, REGRID_NO_MEAN_POSITION

# A field at cell centres, with halos, as saved in the output file; f takes the coordinates of the non-Flat dimensions
# (longitude and latitude in degrees on a LatitudeLongitudeGrid). Built directly from the coordinates, rather than with
# set!, which would compile a kernel for every function.
function centre_data(grid, f)
    underlying = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    horizontal = underlying isa LatitudeLongitudeGrid ? (λnodes, φnodes) : (xnodes, ynodes)
    coordinates = map(n -> isnothing(n) ? [nothing] : vcat(n...), (horizontal[1](underlying, Center(); with_halos = true),
                                                                  horizontal[2](underlying, Center(); with_halos = true),
                                                                  znodes(underlying, Center(); with_halos = true)))
    return [f(filter(!isnothing, (x, y, z))...) for x in coordinates[1], y in coordinates[2], z in coordinates[3]]
end

@testset "Regridding to the mean position" begin

    # Uniform shifts: the label of the grid point y is y - xi, and linear fields are interpolated exactly, so the
    # regridded field is exactly f*(y - xi). A shift of more than half a cell away from a wall leaves the wall row with
    # no mean position (NaN).
    @testset "uniform shift in 3D" begin
        Nx, Ny, Nz = 8, 6, 8
        grid = RectilinearGrid(CPU(), size = (Nx, Ny, Nz), x = (0, 1), y = (0, 1), z = (0, 1),
                               topology = (Bounded, Periodic, Bounded))
        shift = (0.3 / Nx, 0.25 / Ny, 0.7 / Nz) # 0.7 cells up from the bottom
        fstar(x, y, z) = 1 + 2x + 4z
        maps = Tuple(centre_data(grid, (x, y, z) -> s) for s in shift)
        regridded, status, folded = regrid_fields(grid, maps, (centre_data(grid, fstar),))
        xc, yc, zc = xnodes(grid, Center()), ynodes(grid, Center()), znodes(grid, Center())
        exact = [fstar(x - shift[1], y - shift[2], z - shift[3]) for x in xc, y in yc, z in zc]
        @test all(status[:, :, 1] .== REGRID_NO_MEAN_POSITION) && all(isnan, regridded[:, :, 1, 1])
        @test all(status[:, :, 2:end] .== REGRID_OK)
        @test regridded[:, :, 2:end, 1] ≈ exact[:, :, 2:end] atol = 1e-12
        @test folded == 0

        if CUDA.functional()
            gpu_regridded, gpu_status, _ = regrid_fields(grid, maps, (centre_data(grid, fstar),); architecture = GPU())
            @test gpu_status == status
            @test isequal(gpu_regridded[:, :, 1, 1], regridded[:, :, 1, 1]) # the NaNs
            @test gpu_regridded[:, :, 2:end, 1] ≈ regridded[:, :, 2:end, 1] atol = 1e-12
        end
    end

    @testset "uniform shift on a stretched grid" begin
        Nx, Nz = 6, 10
        z_faces = [(k / Nz)^1.5 for k in 0:Nz]
        grid = RectilinearGrid(CPU(), size = (Nx, Nz), x = (0, 1), z = z_faces, topology = (Periodic, Flat, Bounded))
        shift = 0.3 * z_faces[2] # 0.3 of the thinnest cell
        fstar(x, z) = 1 + 4z
        maps = (centre_data(grid, (x, z) -> 0.1 / Nx), nothing, centre_data(grid, (x, z) -> shift))
        regridded, status, _ = regrid_fields(grid, maps, (centre_data(grid, fstar),))
        exact = [fstar(x, z - shift) for x in xnodes(grid, Center()), z in znodes(grid, Center())]
        @test all(status .== REGRID_OK)
        @test regridded[:, 1, :, 1] ≈ exact atol = 1e-12
    end

    # Smooth maps that are zero normal to the walls, against the exact answer (the label found by Newton's method on
    # the analytic map): second order, with no NaNs
    @testset "smooth maps converge at second order" begin
        A = 0.04
        xi(x, z) = (A * sin(2pi * x + 1) * cos(pi * z), A * cos(2pi * x) * sin(pi * z))
        fstar(x, z) = cos(2pi * x) * (1 + z^2) + sin(pi * z)
        function label(y1, y2)
            x, z = y1, y2
            for _ in 1:50
                r1 = x + xi(x, z)[1] - y1; r2 = z + xi(x, z)[2] - y2
                J11 = 1 + 2pi * A * cos(2pi * x + 1) * cos(pi * z); J12 = -pi * A * sin(2pi * x + 1) * sin(pi * z)
                J21 = -2pi * A * sin(2pi * x) * sin(pi * z);       J22 = 1 + pi * A * cos(2pi * x) * cos(pi * z)
                D = J11 * J22 - J12 * J21
                x -= (J22 * r1 - J12 * r2) / D; z -= (-J21 * r1 + J11 * r2) / D
            end
            return x, z
        end
        errors = map((16, 32)) do N
            grid = RectilinearGrid(CPU(), size = (N, N), x = (0, 1), z = (0, 1), topology = (Periodic, Flat, Bounded))
            maps = (centre_data(grid, (x, z) -> xi(x, z)[1]), nothing, centre_data(grid, (x, z) -> xi(x, z)[2]))
            regridded, status, folded = regrid_fields(grid, maps, (centre_data(grid, fstar),))
            @test all(status .== REGRID_OK)
            @test folded == 0
            exact = [fstar(label(x, z)...) for x in xnodes(grid, Center()), z in znodes(grid, Center())]
            maximum(abs, regridded[:, 1, :, 1] .- exact)
        end
        @test errors[2] < errors[1] / 3
    end

    # Compression strong enough that the mean position decreases along x somewhere (1 + a cos(2 pi x) < 0)
    @testset "folding is detected" begin
        a = 1.6
        grid = RectilinearGrid(CPU(), size = (24, 8), x = (0, 1), z = (0, 1), topology = (Periodic, Flat, Bounded))
        maps = (centre_data(grid, (x, z) -> a / 2pi * sin(2pi * x)), nothing, centre_data(grid, (x, z) -> 0))
        _, _, folded = regrid_fields(grid, maps, (centre_data(grid, (x, z) -> x),))
        @test folded > 0
    end

    # A uniform shift in longitude, latitude and depth. The maps are distances in metres, converted to degrees at each
    # cell's latitude, so the eastward map is c_lon degrees of arc along each latitude circle for the shift in longitude
    # to be uniform. A field linear in all three coordinates is then regridded exactly.
    @testset "uniform shift on a LatitudeLongitudeGrid" begin
        Nlon, Nlat, Nz = 8, 6, 4
        grid = LatitudeLongitudeGrid(CPU(), size = (Nlon, Nlat, Nz), longitude = (0, 20), latitude = (30, 50),
                                     z = (-100, 0), topology = (Bounded, Bounded, Bounded))
        R = grid.radius
        c_lon, c_lat, c_z = 0.3 * 20 / Nlon, 0.25 * 20 / Nlat, 0.3 * 100 / Nz # degrees, degrees, metres
        fstar(lon, lat, z) = 1 + 2lon + 3lat + 4z / 100
        maps = (centre_data(grid, (lon, lat, z) -> deg2rad(c_lon) * R * cosd(lat)),
                centre_data(grid, (lon, lat, z) -> deg2rad(c_lat) * R),
                centre_data(grid, (lon, lat, z) -> c_z))
        regridded, status, folded = regrid_fields(grid, maps, (centre_data(grid, fstar),))
        exact = [fstar(lon - c_lon, lat - c_lat, z - c_z) for lon in λnodes(grid, Center()), lat in φnodes(grid, Center()),
                                                             z in znodes(grid, Center())]
        @test all(status .== REGRID_OK)
        @test regridded[:, :, :, 1] ≈ exact atol = 1e-9
        @test folded == 0
    end

    # A flat immersed bottom at z = 0.3: grid points below it are solid, and a shift of 0.7 cells up leaves the lowest
    # fluid row with labels below the bottom (no mean position)
    @testset "immersed bottom" begin
        Nz = 10
        underlying = RectilinearGrid(CPU(), size = (6, Nz), x = (0, 1), z = (0, 1), topology = (Periodic, Flat, Bounded))
        grid = ImmersedBoundaryGrid(underlying, GridFittedBottom(x -> 0.3))
        shift = 0.7 / Nz
        fstar(x, z) = 1 + 4z
        maps = (centre_data(grid, (x, z) -> 0), nothing, centre_data(grid, (x, z) -> shift))
        regridded, status, _ = regrid_fields(grid, maps, (centre_data(grid, fstar),))
        zc = znodes(underlying, Center())
        solid, lowest = findall(<(0.3), zc), findfirst(>(0.3), zc)
        @test all(status[:, 1, solid] .== REGRID_SOLID_POINT)
        @test all(status[:, 1, lowest] .== REGRID_NO_MEAN_POSITION)
        @test all(isnan, regridded[:, 1, 1:lowest, 1])
        exact = [fstar(x, z - shift) for x in xnodes(underlying, Center()), z in zc[lowest+1:end]]
        @test regridded[:, 1, lowest+1:end, 1] ≈ exact atol = 1e-12
    end
end

using NCDatasets
using Oceananigans.Grids: xspacings, yspacings, zspacings, λnodes, φnodes

@testset "jld2_to_netcdf writes interior data, coordinates and spacings" begin
    grids = (
        rectilinear_stretched = RectilinearGrid(size = (6, 4, 5), x = (0, 1), y = (0, 2), z = [-1, -0.6, -0.3, -0.15, -0.05, 0],
                                                topology = (Periodic, Bounded, Bounded)),
        latitude_longitude = LatitudeLongitudeGrid(size = (6, 5, 4), longitude = (0, 60), latitude = (-40, 40), z = (-100, 0)),
        latitude_longitude_stretched = LatitudeLongitudeGrid(size = (4, 5, 3), longitude = [0, 5, 15, 30, 40], latitude = (-40, 40),
                                                             z = [-100, -50, -10, 0]),
    )

    for (label, grid) in pairs(grids)
        @testset "$label" begin
            filename_stem = "test_netcdf_$label"
            model = HydrostaticFreeSurfaceModel(grid; tracers = :c, buoyancy = nothing)
            set!(model, c = (x, y, z) -> x + 10y + 100z, u = (x, y, z) -> y)
            simulation = Simulation(model, Δt = 1, stop_iteration = 2)
            simulation.output_writers[:fields] = JLD2Writer(model, merge(model.velocities, model.tracers);
                                                            filename = filename_stem * ".jld2",
                                                            schedule = IterationInterval(1), overwrite_files = true)
            run!(simulation)

            jld2_to_netcdf(filename_stem * ".jld2", filename_stem * ".nc")

            NCDataset(filename_stem * ".nc") do ds
                # Data: the interior of each saved field, without halos
                for name in ("u", "v", "w", "c")
                    @test ds[name][:, :, :, :] == Array(interior(FieldTimeSeries(filename_stem * ".jld2", name)))
                end

                @test !any(haskey(ds, h) for h in ("Hx", "Hy", "Hz"))

                # Coordinates and spacings: the interior nodes and spacings at each location
                x, y, xnodes_, ynodes_ = grid isa LatitudeLongitudeGrid ? ("lon", "lat", λnodes, φnodes) : ("x", "y", xnodes, ynodes)
                @test ds["$(x)_caa"][:] ≈ xnodes_(grid, Center())
                @test ds["$(x)_faa"][:] ≈ xnodes_(grid, Face())
                @test ds["$(y)_aca"][:] ≈ ynodes_(grid, Center())
                @test ds["$(y)_afa"][:] ≈ ynodes_(grid, Face())
                @test ds["z_aac"][:] ≈ znodes(grid, Center())
                @test ds["z_aaf"][:] ≈ znodes(grid, Face())

                for name in keys(ds)
                    startswith(name, "d") || continue
                    location = Tuple(l == 'f' ? Face() : Center() for l in name[4:end] if l != 'a')
                    spacing = Dict('x' => xspacings, 'y' => yspacings, 'z' => zspacings)[name[2]]
                    @test vec(ds[name][:]) ≈ vec(collect(spacing(grid, location...)))
                end
            end

            rm(filename_stem * ".jld2", force = true)
            rm(filename_stem * ".nc", force = true)
        end
    end
end

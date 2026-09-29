using JLD2
using Oceananigans.Grids: halo_size

@testset "Filtering data saved by an older Oceananigans version" begin
    # reference_sim_oceananigans_v0.110.jld2 is the reference simulation as written by Oceananigans v0.110.
    # Newer versions can't reconstruct the grid (or boundary conditions) saved in it, but the data can be
    # used if the grid of the original simulation is passed to the config.
    old_file = "data/reference_sim_oceananigans_v0.110.jld2"
    kwargs = (; original_data_filename = old_file, var_names_to_filter = ("b",), velocity_names = ("u", "w"),
                N = 2, freq_c = 1e-4 / 4)

    # Without a grid, the grid is read from the file, which fails
    @test_throws Exception OfflineFilterConfig(; kwargs...)

    # The grid of the original simulation, as in reference_sim_offline.jl
    grid = RectilinearGrid(CPU(), size = (10, 10), x = (-5kilometers, 5kilometers), z = (-100, 0),
                           topology = (Periodic, Flat, Bounded))
    config = OfflineFilterConfig(; kwargs..., grid)

    # The data read by the package matches the raw data in the file (which was saved with halos)
    Hx, _, Hz = halo_size(grid)
    b₀, u₀ = jldopen(old_file) do file
        first_iteration = minimum(parse.(Int, filter(!=("serialized"), keys(file["timeseries/b"]))))
        remove_halos(data) = data[Hx+1:end-Hx, :, Hz+1:end-Hz]
        remove_halos(file["timeseries/b/$first_iteration"]), remove_halos(file["timeseries/u/$first_iteration"])
    end

    original_vars = create_original_vars(config)
    @test interior(original_vars.b) == b₀

    reader = create_buffered_reader(config; direction = :forward)
    lo_frame = reader.frames[reader.lo_slot]
    @test interior(lo_frame[:b]) == b₀
    @test interior(lo_frame[:u]) == u₀
end

import CUDA
using Oceananigans.Grids: on_architecture

# ─────────────────────────────────────────────────────────────────────────────
# Shared assertion suite
#
# make_fwd / make_bwd are zero-argument factories that return a freshly
# constructed, primed BufferedDataReader (so each sub-test gets clean state).
# b_frames / u_frames are Vector{Array{Float64,3}} of (nx,1,nz) interior data
# for each stored time, used as ground-truth for interpolation checks.
# ─────────────────────────────────────────────────────────────────────────────
function run_reader_tests(label, make_fwd, make_bwd, model, times, b_frames, u_frames)
    n = length(times)

    @testset "initial buffer state ($label)" begin
        fwd = make_fwd(); bwd = make_bwd()
        @test fwd.direction  == :forward
        @test fwd.T_start    ≈ times[1]
        @test fwd.T_end      ≈ times[n]
        @test fwd.lo_src_idx == 1
        @test fwd.hi_src_idx == 2
        @test bwd.direction  == :backward
        @test bwd.lo_src_idx == n - 1
        @test bwd.hi_src_idx == n
    end

    @testset "forward interpolation at frame boundary ($label)" begin
        fwd = make_fwd()
        # sim_t=0 → physical time = T_start = times[1], α=0, result is frame 1 exactly
        interpolate_to_model!(model, fwd, 0.0)
        @test Array(interior(model.auxiliary_fields.b)) ≈ b_frames[1]
        @test Array(interior(model.velocities.u))       ≈ u_frames[1]
    end

    @testset "forward interpolation at midpoint ($label)" begin
        fwd   = make_fwd()
        t_mid = (times[1] + times[2]) / 2
        sim_t = t_mid - fwd.T_start
        advance_buffer!(fwd, sim_t)
        @test fwd.lo_src_idx == 1
        @test fwd.hi_src_idx == 2
        interpolate_to_model!(model, fwd, sim_t)
        @test Array(interior(model.auxiliary_fields.b)) ≈ 0.5 .* (b_frames[1] .+ b_frames[2])
        @test Array(interior(model.velocities.u))       ≈ 0.5 .* (u_frames[1] .+ u_frames[2])
    end

    @testset "backward interpolation — velocity sign negation ($label)" begin
        bwd   = make_bwd()
        t_mid = (times[n-1] + times[n]) / 2
        sim_t = bwd.T_end - t_mid
        advance_buffer!(bwd, sim_t)
        @test bwd.lo_src_idx == n - 1
        @test bwd.hi_src_idx == n
        interpolate_to_model!(model, bwd, sim_t)
        # Tracer: interpolated, not negated
        @test Array(interior(model.auxiliary_fields.b)) ≈  0.5 .* (b_frames[n-1] .+ b_frames[n])
        # Velocity: interpolated and negated for backward advection
        @test Array(interior(model.velocities.u))       ≈ -0.5 .* (u_frames[n-1] .+ u_frames[n])
    end

    @testset "forward buffer slide ($label)" begin
        fwd = make_fwd()
        @test fwd.lo_src_idx == 1
        @test fwd.hi_src_idx == 2
        # Advance just past times[2]: bracket shifts from [1,2] to [2,3]
        advance_buffer!(fwd, times[2] + 1.0 - fwd.T_start)
        @test fwd.lo_src_idx == 2
        @test fwd.hi_src_idx == 3
    end

    @testset "backward buffer slide ($label)" begin
        bwd = make_bwd()
        @test bwd.lo_src_idx == n - 1
        @test bwd.hi_src_idx == n
        # Physical time just past times[n-2]: bracket shifts from [n-1,n] to [n-2,n-1]
        advance_buffer!(bwd, bwd.T_end - (times[n-2] + 1.0))
        @test bwd.lo_src_idx == n - 2
        @test bwd.hi_src_idx == n - 1
    end

    @testset "advance_buffer! no-op within same bracket ($label)" begin
        fwd = make_fwd()
        lo_slot_before    = fwd.lo_slot
        lo_src_idx_before = fwd.lo_src_idx
        hi_src_idx_before = fwd.hi_src_idx
        # Still inside [times[1], times[2]] — state must be unchanged
        advance_buffer!(fwd, times[2] - 1.0 - fwd.T_start)
        @test fwd.lo_slot    == lo_slot_before
        @test fwd.lo_src_idx == lo_src_idx_before
        @test fwd.hi_src_idx == hi_src_idx_before
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# JLD2 source
# ─────────────────────────────────────────────────────────────────────────────
@testset "BufferedDataReader — JLD2 source" begin

    datafile = "data/reference_sim.jld2"
    b_fts    = FieldTimeSeries(datafile, "b")
    u_fts    = FieldTimeSeries(datafile, "u")
    grid     = b_fts.grid
    times    = Float64.(b_fts.times)
    n        = length(times)

    has_gpu    = try; on_architecture(GPU(), grid); true; catch; false; end
    test_archs = has_gpu ? [CPU(), GPU()] : [CPU()]

    for arch in test_archs
        arch_label = arch isa CPU ? "CPU" : "GPU"
        arch_grid  = on_architecture(arch, grid)

        config = OfflineFilterConfig(
            original_data_filename  = datafile,
            var_names_to_filter     = ("b",),
            velocity_names          = ("u", "w"),
            N = 1, freq_c = 1e-4,
            architecture = arch,
            grid = grid,
            map_to_mean             = false,
            output_netcdf           = false,
            compute_Eulerian_filter = false,
            delete_intermediate_files = false,
        )

        model = NonhydrostaticModel(arch_grid;
                    auxiliary_fields = (b = Field{Center,Center,Center}(arch_grid),))

        # Reference data as plain CPU arrays
        b_frames = [Array(interior(b_fts[i])) for i in 1:n]
        u_frames = [Array(interior(u_fts[i])) for i in 1:n]

        make_fwd = () -> create_buffered_reader(config; direction = :forward)
        make_bwd = () -> create_buffered_reader(config; direction = :backward)

        run_reader_tests(arch_label, make_fwd, make_bwd, model, times, b_frames, u_frames)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# NetCDF source
# ─────────────────────────────────────────────────────────────────────────────
@testset "BufferedDataReader — NetCDF source" begin

    datafile = "data/reference_sim.nc"
    b_fts    = FieldTimeSeries(datafile, "b")
    u_fts    = FieldTimeSeries(datafile, "u")
    grid     = b_fts.grid
    times    = Float64.(b_fts.times)
    n        = length(times)

    has_gpu    = try; on_architecture(GPU(), grid); true; catch; false; end
    test_archs = has_gpu ? [CPU(), GPU()] : [CPU()]

    # ── NetCDFDataSource-specific tests ──────────────────────────────────────

    @testset "NetCDFDataSource construction" begin
        src = NetCDFDataSource(datafile, ["b"], ["u"])
        @test src.stored_times   ≈ times
        @test src.stored_indices == collect(1:n)

        # T_start / T_end windowing selects a sub-range of frames
        src_w = NetCDFDataSource(datafile, ["b"], ["u"]; T_start = 7200.0, T_end = 14400.0)
        @test src_w.stored_times   ≈ [7200.0, 10800.0, 14400.0]
        @test src_w.stored_indices == [3, 4, 5]
    end

    @testset "read_frame!" begin
        src   = NetCDFDataSource(datafile, ["b"], ["u"])
        b_cpu = Field{Center,Center,Center}(grid)

        OceananigansLagrangianFilter.DataIO.read_frame!(b_cpu, src, "b", 1)
        @test interior(b_cpu) ≈ interior(b_fts[1])

        OceananigansLagrangianFilter.DataIO.read_frame!(b_cpu, src, "b", 4)
        @test interior(b_cpu) ≈ interior(b_fts[4])
    end

    # ── Shared reader-behaviour tests (same as JLD2) ──────────────────────────

    for arch in test_archs
        arch_label = arch isa CPU ? "CPU" : "GPU"
        arch_grid  = on_architecture(arch, grid)

        config = OfflineFilterConfig(
            original_data_filename  = datafile,
            var_names_to_filter     = ("b",),
            velocity_names          = ("u", "w"),
            N = 1, freq_c = 1e-4,
            architecture = arch,
            grid = grid,
            map_to_mean             = false,
            output_netcdf           = false,
            compute_Eulerian_filter = false,
            delete_intermediate_files = false,
        )

        model = NonhydrostaticModel(arch_grid;
                    auxiliary_fields = (b = Field{Center,Center,Center}(arch_grid),))

        b_frames = [Array(interior(b_fts[i])) for i in 1:n]
        u_frames = [Array(interior(u_fts[i])) for i in 1:n]

        make_fwd = () -> create_buffered_reader(config; direction = :forward)
        make_bwd = () -> create_buffered_reader(config; direction = :backward)

        run_reader_tests(arch_label, make_fwd, make_bwd, model, times, b_frames, u_frames)
    end
end

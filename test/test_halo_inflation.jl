using Oceananigans.Grids: halo_size

@testset "shared_halo_regions with inflated model halo" begin
    # shared_halo_regions is internal (not exported). Copy through the regions it returns, as interpolate_to_model! does.
    shared_halo_regions = OceananigansLagrangianFilter.DataIO.shared_halo_regions
    function copy_input_data!(field, data_field)
        field_region, data_region = shared_halo_regions(field, data_field)
        view(parent(field), field_region...) .= view(parent(data_field), data_region...)
    end

    Nx, Nz = 6, 5
    topology = (Periodic, Flat, Bounded)
    data_grid  = RectilinearGrid(size = (Nx, Nz), x = (0, 1), z = (-1, 0), halo = (3, 3); topology)
    model_grid = RectilinearGrid(size = (Nx, Nz), x = (0, 1), z = (-1, 0), halo = (5, 5); topology)

    # Test each location used by the offline filter: u, w, and tracers/auxiliary fields
    for loc in ((Face(), Center(), Center()), (Center(), Center(), Face()), (Center(), Center(), Center()))
        data_field = Field(loc, data_grid)
        # Distinct, non-zero value in every cell, including halos
        parent(data_field) .= reshape(1:length(parent(data_field)), size(parent(data_field)))

        # Matching halos: identical to copying the whole parent array
        same_halo_field = Field(loc, data_grid)
        copy_input_data!(same_halo_field, data_field)
        @test parent(same_halo_field) == parent(data_field)

        # Inflated model halo: the interior and the 3 shared halo cells are copied...
        model_field = Field(loc, model_grid)
        copy_input_data!(model_field, data_field)
        @test interior(model_field) == interior(data_field)
        @test model_field.data[axes(data_field.data)...] == data_field.data
        # ...and the 2 extra outer halo layers are left as zero
        @test count(!iszero, parent(model_field)) == length(parent(data_field))
    end
end

@testset "Offline run test: WENO9 with data saved at a smaller halo" begin
    # reference_sim.jld2 was saved from a WENO5 simulation, so WENO(order=9) makes the LagrangianFilter
    # inflate its grid halo in both a Periodic (x) and a Bounded (z) direction. This used to throw a
    # DimensionMismatch when copying the saved data into the model fields.
    # There's no reference output for WENO9, so we only check the run completes and gives finite, non-trivial output.
    # Kept short (T = 3hours, no maps) since most of the cost is compiling the WENO9 kernels.
    data_halo = halo_size(FieldTimeSeries("data/reference_sim.jld2", "u").grid)
    @test data_halo[1] < 5 && data_halo[3] < 5 # Make sure this test actually exercises halo inflation

    test_filename_stem = "data/test_offline_output_WENO9"
    try
        filter_config = OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                        output_filename = test_filename_stem * ".jld2",
                                        var_names_to_filter = ("b",),
                                        velocity_names = ("u", "w"),
                                        architecture = CPU(),
                                        Δt = 20minutes,
                                        T = 3hours,
                                        T_out = 1hour,
                                        N = 2,
                                        freq_c = 1e-4 / 4,
                                        advection = WENO(order = 9), # Requires halo 5, larger than the saved data
                                        compute_mean_velocities = false,
                                        map_to_mean = false,
                                        delete_intermediate_files = true)

        run_offline_Lagrangian_filter(filter_config)

        @test isfile(test_filename_stem * ".jld2")
        b_filtered = FieldTimeSeries(test_filename_stem * ".jld2", "b_Lagrangian_filtered")
        @test all(isfinite, b_filtered.data)
        @test any(!=(0), b_filtered.data)
    finally
        rm(test_filename_stem * ".jld2", force = true)
    end
end

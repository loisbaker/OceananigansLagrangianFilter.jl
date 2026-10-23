# Shared helpers used across the integration tests in test_offline_integration.jl
# and test_online_integration.jl.

using JLD2
using NCDatasets

# The architectures to run the integration tests against reference data on: the GPU too, if CUDA finds one
test_architectures() = CUDA.functional() ? (CPU(), GPU()) : (CPU(),)

# "CPU" or "GPU", for test names and output filenames
arch_label(arch) = string(nameof(typeof(arch)))

"""
    compare_filter_output_to_reference(test_filename_stem, ref_filename, varname; rtol = 1e-6)

Compares `varname` and its `_Lagrangian_filtered`, `_Lagrangian_filtered_at_mean`, and
`_Eulerian_filtered` variants between a freshly generated `<test_filename_stem>.jld2` and
a saved reference file, then deletes the generated `.jld2`/`.nc` files.

The explicit `rtol` (rather than relying on `isapprox`'s default tolerance) is so that
small floating-point differences introduced by, e.g., a new Julia/Oceananigans/BLAS
version don't silently change what "close enough" means between CI runs.
"""
function compare_filter_output_to_reference(test_filename_stem::String, ref_filename::String,
                                            varname::String; rtol = 1e-6)
    test_filename = test_filename_stem * ".jld2"

    @test isfile(test_filename)
    @test isfile(test_filename_stem * ".nc")

    for suffix in ("", "_Lagrangian_filtered", "_Lagrangian_filtered_at_mean", "_Eulerian_filtered")
        name = varname * suffix
        test_fts = FieldTimeSeries(test_filename, name)
        ref_fts = FieldTimeSeries(ref_filename, name)
        @test isapprox(test_fts.data, ref_fts.data; rtol)
        @test test_fts.grid == ref_fts.grid

        # The NetCDF copy holds the interior of each field only, without halo cells
        NCDataset(test_filename_stem * ".nc") do ds
            @test ds[name][:, :, :, :] == Array(interior(test_fts))
        end
    end

    return nothing
end

"""
    reference_output_without_Eulerian_filter(filename)

Return the path to a temporary copy of the filter output `filename` without its Eulerian filtered outputs, for testing
`compute_Eulerian_filter!` on it. The original file is not changed.
"""
function reference_output_without_Eulerian_filter(filename)
    filename_copy = joinpath(mktempdir(), basename(filename))
    cp(filename, filename_copy)
    jldopen(filename_copy, "a+") do file
        for name in filter(contains("Eulerian"), collect(keys(file["timeseries"])))
            Base.delete!(file["timeseries"], name)
        end
    end
    return filename_copy
end

"""
    reference_sim_with_extra_variable()

Return the path to a temporary copy of `data/reference_sim.jld2` with an extra variable `b2 = 2b`, for testing
`extra_auxiliary_var_names`. The reference data itself is not changed.
"""
function reference_sim_with_extra_variable()
    filename = joinpath(mktempdir(), "reference_sim_extra.jld2")
    cp("data/reference_sim.jld2", filename)
    jldopen(filename, "a+") do file
        b = file["timeseries/b"]
        for key in keys(b["serialized"]) # location, indices and boundary conditions
            file["timeseries/b2/serialized/$key"] = b["serialized/$key"]
        end
        for key in filter(!=("serialized"), keys(b)) # one array per saved iteration
            file["timeseries/b2/$key"] = 2 .* b[key]
        end
    end
    return filename
end

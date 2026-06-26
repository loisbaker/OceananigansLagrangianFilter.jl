# Shared helpers used across the integration tests in test_offline_integration.jl
# and test_online_integration.jl.

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
    end

    return nothing
end

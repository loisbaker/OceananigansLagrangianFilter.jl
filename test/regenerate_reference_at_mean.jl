# Regenerate only the regridded variables (with the suffix _at_mean) in the reference outputs, after a change to
# regrid_to_mean_position!, without rerunning the reference simulations. Everything else in each reference file is
# copied unchanged, and checked to be unchanged, before the file is replaced. Run from the test directory, after
# reference_sim_offline.jl and reference_sim_online.jl have made the reference outputs:
#
#     julia --project=.. regenerate_reference_at_mean.jl

using OceananigansLagrangianFilter
using Oceananigans
using Oceananigans.Units
using JLD2
using JLD2: Group
using Printf

f = 1e-4 # Coriolis frequency [s^-1], as in the reference simulations

# The filter configurations of reference_sim_offline.jl and reference_sim_online.jl, writing to `filename`
offline_config(filename) = OfflineFilterConfig(original_data_filename = "data/reference_sim.jld2",
                                               output_filename = filename,
                                               var_names_to_filter = ("b",),
                                               velocity_names = ("u", "w"),
                                               architecture = CPU(),
                                               Δt = 20minutes,
                                               T_out = 1hour,
                                               N = 2,
                                               freq_c = f/4,
                                               compute_mean_velocities = true,
                                               output_netcdf = false,
                                               delete_intermediate_files = true,
                                               compute_Eulerian_filter = true)

online_config(filename) = OnlineFilterConfig(grid = jldopen(file -> file["serialized/grid"], filename),
                                             output_filename = filename,
                                             var_names_to_filter = ("b", "T"),
                                             velocity_names = ("u", "w"),
                                             N = 1,
                                             freq_c = f/4)

# Copy the group `source` into `destination`, leaving out the regridded variables. Writing them into a fresh file,
# rather than replacing them in place, stops the file growing (JLD2 does not reuse the space of deleted data).
function copy_without_at_mean!(destination, source)
    for key in keys(source)
        endswith(key, "_at_mean") && continue
        object = source[key]
        if object isa JLD2.Group
            copy_without_at_mean!(Group(destination, key), object)
        else
            destination[key] = object
        end
    end
end

# The paths of every dataset in the group, except the regridded variables
function dataset_paths(group, prefix = "")
    paths = String[]
    for key in keys(group)
        endswith(key, "_at_mean") && continue
        path = prefix * key
        group[key] isa JLD2.Group ? append!(paths, dataset_paths(group[key], path * "/")) : push!(paths, path)
    end
    return paths
end

same(a, b) = isequal(a, b) || repr(a) == repr(b)

for (reference, make_config) in (("data/reference_offline_output.jld2", offline_config),
                                 ("data/reference_online_output.jld2", online_config))
    new_file = joinpath(mktempdir(), basename(reference))
    jldopen(source -> jldopen(destination -> copy_without_at_mean!(destination, source), new_file, "w"), reference)
    regrid_to_mean_position!(make_config(new_file))

    jldopen(reference) do old
        jldopen(new_file) do new
            # Everything but the regridded variables is unchanged
            paths = dataset_paths(old)
            changed = filter(path -> !haskey(new, path) || !same(old[path], new[path]), paths)
            isempty(changed) || error("$reference: these would change: $changed")
            @info @sprintf("%s: %d datasets copied unchanged", reference, length(paths))

            # The regridded variables have no NaNs (the reference cases have no immersed boundary and flat walls), and
            # differ from the previous ones by a fraction of the field's variation
            for name in filter(endswith("_at_mean"), collect(keys(old["timeseries"])))
                iterations = filter(!=("serialized"), collect(keys(old["timeseries/$name"])))
                nans = sum(count(isnan, new["timeseries/$name/$i"]) for i in iterations)
                nans == 0 || error("$reference: $name has $nans NaNs")
                filtered = replace(name, "_at_mean" => "")
                difference = maximum(iterations) do i
                    field = old["timeseries/$filtered/$i"]
                    variation = maximum(abs, field .- sum(field) / length(field))
                    variation == 0 ? 0.0 : maximum(abs, new["timeseries/$name/$i"] .- old["timeseries/$name/$i"]) / variation
                end
                @info @sprintf("   %s: no NaNs; largest change %.1e of the field's variation", name, difference)
            end
        end
    end

    @info @sprintf("%s: %d bytes before, %d after", reference, filesize(reference), filesize(new_file))
    mv(new_file, reference; force = true)
end

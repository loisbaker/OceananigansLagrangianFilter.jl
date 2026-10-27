# Running an offline filter in batches of its terms (the batch_size option of run_offline_Lagrangian_filter): each batch
# runs the forward and backward passes with a smaller model, and its output is merged into the output of the whole filter

# ── Batches of filter terms ───────────────────────────────────────────────────

# Run the filter in `config` in batches of `batch_size` terms, merging each batch's output into `config.output_filename`
function run_batches(config, batch_size)
    batch_size >= 1 || error("batch_size must be at least 1, got $batch_size")
    N_terms = config.filter_params.N_coeffs
    batches = [i:min(i + batch_size - 1, N_terms) for i in 1:batch_size:N_terms]

    # Each batch builds its model and outputs from the same state of the counter that Oceananigans uses to name
    # interpolation operators (see run_forward_and_backward_passes), so that batches of the same size have the same
    # types and are only compiled once
    identity_counter = output_identity_counter()

    for (n, terms) in enumerate(batches)
        @info "Running batch $n of $(length(batches)): filter terms $(first(terms)) to $(last(terms))"
        batch_config = batch_filter_config(config, terms, n)
        restore_output_identity_counter!(identity_counter)
        run_forward_and_backward_passes(batch_config)

        merge_batch_output!(config.output_filename, batch_config, terms, n)
        config.delete_intermediate_files && rm(batch_config.output_filename)
        @info "Merged batch $n into $(config.output_filename)"
    end

    return nothing
end

# The config for the batch of filter terms `terms`, numbered from 1, writing to its own files, without post-processing
function batch_filter_config(config, terms, n)
    filter_params = config.filter_params
    batch_coefficients = NamedTuple(Symbol(x, j) => getproperty(filter_params, Symbol(x, i))
                                    for (j, i) in enumerate(terms) for x in ("a", "b", "c", "d"))
    other_params = NamedTuple(k => v for (k, v) in pairs(filter_params) if !occursin(r"^[abcd][0-9]+$", String(k)))
    batch_filter_params = merge(other_params, batch_coefficients, (; N_coeffs = length(terms)))

    batch_filename(filename) = filename[1:end-5] * "_batch$n.jld2"
    replacements = (; filter_params = batch_filter_params,
                      output_filename = batch_filename(config.output_filename),
                      forward_output_filename = batch_filename(config.forward_output_filename),
                      backward_output_filename = batch_filename(config.backward_output_filename),
                      regrid_to_mean = false, compute_Eulerian_filter = false, output_netcdf = false)

    return OfflineFilterConfig((get(replacements, name, getfield(config, name)) for name in fieldnames(OfflineFilterConfig))...)
end

# Merge the output of batch `n` of filter terms `terms` into `output_filename`. Separate outputs of each term
# (<name>_C<j>_scaled and <name>_S<j>_scaled) are renumbered to the term's number in the whole filter. Combined outputs
# are added to the sum of the batches so far, written to a new file that then replaces `output_filename` (as JLD2
# doesn't reclaim the space of overwritten data). The other outputs (the time and original data) and the file
# metadata are the same in every batch, so they are copied from the first
function merge_batch_output!(output_filename, batch_config, terms, n)
    outputs, _ = filter_output_options(batch_config.filter_params)
    filtered_names = filtered_output_names(batch_config)

    if outputs === :separate
        jldopen(output_filename, n == 1 ? "w" : "r+") do merged_file
            jldopen(batch_config.output_filename, "r") do batch_file
                n == 1 && copy_batch_shared_outputs!(merged_file, batch_file, filtered_names)
                for name in filtered_names
                    stem, kind, j = match(r"^(.*)_([CS])([0-9]+)_scaled$", name).captures
                    merged_name = "$(stem)_$(kind)$(terms[parse(Int, j)])_scaled"
                    copy_jld2_path!(batch_file, "timeseries/$name", merged_file, "timeseries/$merged_name")
                end
            end
        end
    elseif n == 1
        cp(batch_config.output_filename, output_filename; force = true)
    else
        summed_filename = output_filename * ".merging.jld2"
        jldopen(summed_filename, "w") do summed_file
            jldopen(output_filename, "r") do merged_file
                jldopen(batch_config.output_filename, "r") do batch_file
                    copy_batch_shared_outputs!(summed_file, merged_file, filtered_names)
                    for name in filtered_names
                        copy_jld2_path!(merged_file, "timeseries/$name/serialized", summed_file, "timeseries/$name/serialized")
                        for iter in keys(merged_file["timeseries/$name"])
                            iter == "serialized" && continue
                            summed_file["timeseries/$name/$iter"] = merged_file["timeseries/$name/$iter"] .+
                                                                    batch_file["timeseries/$name/$iter"]
                        end
                    end
                end
            end
        end
        mv(summed_filename, output_filename; force = true)
    end

    return nothing
end

# Copy the outputs shared by all the batches (the file metadata, the time and the original data), i.e. everything in
# `source` apart from the filtered outputs `filtered_names`, to `dest`
function copy_batch_shared_outputs!(dest, source, filtered_names)
    for key in keys(source)
        key == "timeseries" || copy_jld2_path!(source, key, dest, key)
    end
    for name in keys(source["timeseries"])
        name in filtered_names || copy_jld2_path!(source, "timeseries/$name", dest, "timeseries/$name")
    end
    return nothing
end

# Copy the dataset or group at `source_path` in `source` to `dest_path` in `dest`
function copy_jld2_path!(source, source_path, dest, dest_path)
    if source[source_path] isa JLD2.Group
        JLD2.Group(dest, dest_path)
        for key in keys(source[source_path])
            copy_jld2_path!(source, "$source_path/$key", dest, "$dest_path/$key")
        end
    else
        dest[dest_path] = source[source_path]
    end
    return nothing
end

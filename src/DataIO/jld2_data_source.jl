"""
    JLD2DataSource

Wraps a JLD2 file of Oceananigans output, providing indexed frame reads without
writing any intermediate files.

Fields
======
- `filename`           : path to the JLD2 file
- `var_names`          : tracer/auxiliary variable names to read
- `vel_names`          : velocity variable names to read
- `stored_times`       : physical times for the T_start..T_end window (ascending)
- `stored_iterations`  : corresponding JLD2 iteration keys
"""
struct JLD2DataSource <: AbstractDataSource
    filename         :: String
    var_names        :: Vector{String}
    vel_names        :: Vector{String}
    stored_times     :: Vector{Float64}
    stored_iterations:: Vector{Int}
end

"""
    JLD2DataSource(filename, var_names, vel_names; T_start=nothing, T_end=nothing)

Open `filename` and build a source covering the time window [`T_start`, `T_end`].
"""
function JLD2DataSource(filename::String,
                        var_names,
                        vel_names;
                        T_start = nothing,
                        T_end   = nothing)

    iters, times = jldopen(filename, "r") do f
        raw_iters = sort!(parse.(Int, keys(f["timeseries/t"])))
        raw_times = Float64[f["timeseries/t/$i"] for i in raw_iters]
        raw_iters, raw_times
    end

    lo = isnothing(T_start) ? 1              : searchsortedfirst(times, T_start)
    hi = isnothing(T_end)   ? length(times)  : searchsortedlast(times, T_end)

    return JLD2DataSource(filename,
                          collect(String, var_names),
                          collect(String, vel_names),
                          times[lo:hi],
                          iters[lo:hi])
end

"""
    read_frame!(cpu_field, source::JLD2DataSource, varname, idx)

Read frame `idx` of `varname` from the JLD2 file into `cpu_field` (a CPU-resident
Oceananigans `Field`). Fills the interior, then fills halo regions.
"""
function read_frame!(cpu_field::Field, source::JLD2DataSource, varname::String, idx::Int)
    iter = source.stored_iterations[idx]
    data = jldopen(source.filename, "r") do f
        Array(f["timeseries/$varname/$iter"])
    end

    # The JLD2 writer saves interior-only data by default. Handle the case where
    # the file was written with halos (parent-size data) for robustness.
    f_interior = interior(cpu_field)
    if size(data) == size(parent(cpu_field))
        parent(cpu_field) .= data
    else
        f_interior .= reshape(data, size(f_interior))
        fill_halo_regions!(cpu_field)
    end
    return nothing
end

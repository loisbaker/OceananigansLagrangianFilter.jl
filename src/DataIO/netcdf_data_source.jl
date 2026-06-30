"""
    NetCDFDataSource

Wraps a NetCDF file of Oceananigans (or compatible) output, providing indexed
frame reads without writing any intermediate files.

Fields
======
- `filename`       : path to the NetCDF file
- `var_names`      : tracer/auxiliary variable names to read
- `vel_names`      : velocity variable names to read
- `stored_times`   : physical times for the T_start..T_end window (ascending)
- `stored_indices` : corresponding indices into the NetCDF time dimension
- `time_name`      : name of the time coordinate in the file (default `"time"`)
"""
struct NetCDFDataSource <: AbstractDataSource
    filename      :: String
    var_names     :: Vector{String}
    vel_names     :: Vector{String}
    stored_times  :: Vector{Float64}
    stored_indices:: Vector{Int}
    time_name     :: String
end

"""
    NetCDFDataSource(filename, var_names, vel_names;
                     T_start=nothing, T_end=nothing, time_name="time")

Open `filename` and build a source covering the time window [`T_start`, `T_end`].
"""
function NetCDFDataSource(filename::String,
                          var_names,
                          vel_names;
                          T_start   = nothing,
                          T_end     = nothing,
                          time_name = "time")

    times, file_indices = NCDatasets.Dataset(filename, "r") do ds
        raw_times = Float64.(ds[time_name][:])
        collect(Float64, raw_times), collect(Int, 1:length(raw_times))
    end

    lo = isnothing(T_start) ? 1              : searchsortedfirst(times, T_start)
    hi = isnothing(T_end)   ? length(times)  : searchsortedlast(times, T_end)

    return NetCDFDataSource(filename,
                            collect(String, var_names),
                            collect(String, vel_names),
                            times[lo:hi],
                            file_indices[lo:hi],
                            time_name)
end

"""
    read_frame!(cpu_field, source::NetCDFDataSource, varname, idx)

Read frame `idx` (in the source's time window) of `varname` from the NetCDF file
into `cpu_field`. Fills the interior, then fills halo regions.

Assumes Oceananigans NetCDF output where the time dimension is last.
"""
function read_frame!(cpu_field::Field, source::NetCDFDataSource, varname::String, idx::Int)
    file_idx = source.stored_indices[idx]

    data = NCDatasets.Dataset(source.filename, "r") do ds
        v = ds[varname]
        ndim = ndims(v) - 1   # all spatial dims; time is last
        # Build index: (:, :, ..., file_idx)
        idx_tuple = (ntuple(_ -> Colon(), ndim)..., file_idx)
        Array(v[idx_tuple...])
    end

    f_interior = interior(cpu_field)
    f_interior .= reshape(data, size(f_interior))
    fill_halo_regions!(cpu_field)
    return nothing
end

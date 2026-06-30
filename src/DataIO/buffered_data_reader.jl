"""
    BufferedDataReader

Reads simulation output frame-by-frame from disk into a two-slot GPU buffer,
interpolating in time on demand. Replaces the intermediate `_filter_input.jld2`
file and the `FieldTimeSeries`-based input pipeline.

Two buffer slots each hold a full set of GPU `Field`s (one per variable). As
the simulation advances, the slot holding the older frame is recycled: the swap
is just an integer flip (`lo_slot = 3 - lo_slot`) with zero allocation.

Forward direction  : physical time = T_start + sim_t  (buffer scans file forward)
Backward direction : physical time = T_end   - sim_t  (buffer scans file backward;
                     velocity signs are negated during interpolation)
"""
mutable struct BufferedDataReader{S <: AbstractDataSource, F, C}
    source    :: S
    direction :: Symbol          # :forward or :backward
    T_start   :: Float64
    T_end     :: Float64
    vel_names :: Vector{String}
    var_names :: Vector{String}

    # frames[slot][varname] → GPU Field
    frames     :: F
    # cpu_frames[slot][varname] → CPU Field (staging area for disk reads)
    cpu_frames :: C

    lo_slot    :: Int   # which slot (1 or 2) holds the lo-time frame
    lo_src_idx :: Int   # index into source.stored_times for the lo frame
    hi_src_idx :: Int   # index into source.stored_times for the hi frame
end

# ──────────────────────────────────────────────────────────────────────────────
# Time mapping helpers
# ──────────────────────────────────────────────────────────────────────────────

physical_time(r::BufferedDataReader, sim_t) =
    r.direction == :forward ? r.T_start + sim_t : r.T_end - sim_t

function bracket_indices(r::BufferedDataReader, sim_t)
    t_phys = physical_time(r, sim_t)
    times  = stored_times(r.source)
    # Find lo s.t. times[lo] ≤ t_phys ≤ times[lo+1]
    lo = clamp(searchsortedlast(times, t_phys), 1, length(times) - 1)
    return lo, lo + 1
end

function _interp_weight(r::BufferedDataReader, sim_t)
    t_phys = physical_time(r, sim_t)
    times  = stored_times(r.source)
    t_lo   = times[r.lo_src_idx]
    t_hi   = times[r.hi_src_idx]
    return (t_phys - t_lo) / (t_hi - t_lo)
end

# ──────────────────────────────────────────────────────────────────────────────
# Frame loading
# ──────────────────────────────────────────────────────────────────────────────

function _load_frame!(r::BufferedDataReader, slot::Int, src_idx::Int)
    all_names = [r.vel_names; r.var_names]
    for varname in all_names
        vsym      = Symbol(varname)
        cpu_field = r.cpu_frames[slot][vsym]
        gpu_field = r.frames[slot][vsym]

        read_frame!(cpu_field, r.source, varname, src_idx)

        # Bulk host → device transfer (single copy of parent including halos)
        copyto!(parent(gpu_field), parent(cpu_field))
    end
    return nothing
end

# ──────────────────────────────────────────────────────────────────────────────
# Buffer advancement (sliding-window with swap)
# ──────────────────────────────────────────────────────────────────────────────

"""
    advance_buffer!(reader, sim_t)

Ensure the two buffer slots bracket `sim_t`. Reuses whichever frame is still
valid (forward or backward slide) so at most one new frame is read per call.
"""
function advance_buffer!(r::BufferedDataReader, sim_t)
    lo_new, hi_new = bracket_indices(r, sim_t)
    (lo_new == r.lo_src_idx && hi_new == r.hi_src_idx) && return

    if lo_new == r.hi_src_idx
        # ── Forward slide: old hi becomes new lo ─────────────────────────────
        r.lo_slot    = 3 - r.lo_slot    # swap: old hi slot is now lo slot
        r.lo_src_idx = lo_new
        r.hi_src_idx = hi_new
        _load_frame!(r, 3 - r.lo_slot, hi_new)   # load new hi into free slot

    elseif hi_new == r.lo_src_idx
        # ── Backward slide: old lo becomes new hi ────────────────────────────
        new_lo_slot  = 3 - r.lo_slot   # old hi slot will hold new lo
        _load_frame!(r, new_lo_slot, lo_new)
        r.lo_slot    = new_lo_slot
        r.lo_src_idx = lo_new
        r.hi_src_idx = hi_new

    else
        # ── Jump (initialisation or large time skip) ─────────────────────────
        r.lo_src_idx = lo_new
        r.hi_src_idx = hi_new
        _load_frame!(r, r.lo_slot,         lo_new)
        _load_frame!(r, 3 - r.lo_slot,     hi_new)
    end
    return nothing
end

# ──────────────────────────────────────────────────────────────────────────────
# Interpolation into model fields
# ──────────────────────────────────────────────────────────────────────────────

"""
    interpolate_to_model!(model, reader, sim_t)

Linearly interpolate the two buffered frames at `sim_t` and write the result
directly into `model.velocities` and `model.auxiliary_fields`, copying the
full parent array (interior + halos) as the existing pipeline does.

Velocities are negated for backward-direction filtering.
"""
function interpolate_to_model!(model, r::BufferedDataReader, sim_t)
    α        = _interp_weight(r, sim_t)
    β        = 1 - α
    vel_sign = r.direction == :backward ? -1 : 1

    lo = r.frames[r.lo_slot]
    hi = r.frames[3 - r.lo_slot]

    for vname in r.vel_names
        vsym = Symbol(vname)
        dest = getproperty(model.velocities, vsym)
        parent(dest) .= vel_sign .* (β .* parent(lo[vsym]) .+ α .* parent(hi[vsym]))
    end

    for vname in r.var_names
        vsym = Symbol(vname)
        dest = getproperty(model.auxiliary_fields, vsym)
        parent(dest) .= β .* parent(lo[vsym]) .+ α .* parent(hi[vsym])
    end

    return nothing
end

# ──────────────────────────────────────────────────────────────────────────────
# Buffer field allocation helpers
# ──────────────────────────────────────────────────────────────────────────────

# Build one CPU template and one GPU template per variable, using a single
# FieldTimeSeries call per variable (JLD2) or the config grid (NetCDF).
# Returns (cpu_templates, gpu_templates) as Dict{Symbol, Field}.
function _build_templates(source::JLD2DataSource, all_names, arch)
    cpu_templates = Dict{Symbol, Field}()
    gpu_templates = Dict{Symbol, Field}()
    for v in all_names
        vsym = Symbol(v)
        # One FieldTimeSeries per variable — reads metadata only, no frame data loaded.
        fts = FieldTimeSeries(source.filename, v; architecture = CPU(), backend = InMemory(2))
        loc = Oceananigans.Fields.location(fts)
        cpu_templates[vsym] = Field{loc[1], loc[2], loc[3]}(fts.grid)
        gpu_templates[vsym] = Field{loc[1], loc[2], loc[3]}(on_architecture(arch, fts.grid))
    end
    return cpu_templates, gpu_templates
end

function _build_templates(source::NetCDFDataSource, all_names, arch, grid, locations)
    cpu_templates = Dict{Symbol, Field}()
    gpu_templates = Dict{Symbol, Field}()
    for v in all_names
        vsym = Symbol(v)
        if haskey(locations, v)
            loc = locations[v]
            cpu_templates[vsym] = Field{loc[1], loc[2], loc[3]}(on_architecture(CPU(), grid))
            gpu_templates[vsym] = Field{loc[1], loc[2], loc[3]}(on_architecture(arch,  grid))
        else
            # Infer location and grid from the NetCDF file attributes (written by Oceananigans NetCDFWriter)
            fts = FieldTimeSeries(source.filename, v; architecture = CPU(), backend = InMemory(2))
            loc = Oceananigans.Fields.location(fts)
            cpu_templates[vsym] = Field{loc[1], loc[2], loc[3]}(fts.grid)
            gpu_templates[vsym] = Field{loc[1], loc[2], loc[3]}(on_architecture(arch, fts.grid))
        end
    end
    return cpu_templates, gpu_templates
end

# ──────────────────────────────────────────────────────────────────────────────
# Public constructor
# ──────────────────────────────────────────────────────────────────────────────

"""
    create_buffered_reader(config; direction=:forward, locations=Dict(), time_name="time")

Build a `BufferedDataReader` for the data file named in `config`, covering the
time window [`config.T_start`, `config.T_end`].

Keyword arguments
=================
- `direction`  : `:forward` (default) or `:backward`
- `locations`  : `Dict{String, NTuple{3, DataType}}` mapping variable names to their
                 spatial location, e.g. `Dict("u" => (Face,Center,Center))`.
                 Only required for NetCDF sources (JLD2 sources infer location from
                 file metadata). Defaults to `(Center,Center,Center)` for all vars.
- `time_name`  : name of the time coordinate in a NetCDF file (default `"time"`)
"""
function create_buffered_reader(config::AbstractConfig;
                                direction :: Symbol                           = :forward,
                                locations :: Dict{String, NTuple{3, DataType}} = Dict{String, NTuple{3, DataType}}(),
                                time_name :: String                            = "time")

    filename  = config.original_data_filename
    var_names = collect(String, config.var_names_to_filter)
    vel_names = collect(String, config.velocity_names)
    arch      = config.architecture
    all_names = [vel_names; var_names]

    # Build the data source ────────────────────────────────────────────────────
    source = if endswith(filename, ".nc") || endswith(filename, ".netcdf")
        NetCDFDataSource(filename, var_names, vel_names;
                         T_start   = Float64(config.T_start),
                         T_end     = Float64(config.T_end),
                         time_name = time_name)
    else
        JLD2DataSource(filename, var_names, vel_names;
                       T_start = Float64(config.T_start),
                       T_end   = Float64(config.T_end))
    end

    if length(stored_times(source)) < 2
        error("BufferedDataReader requires at least 2 frames in [T_start, T_end]. " *
              "Found $(length(stored_times(source))).")
    end

    # Build one template per variable (one FieldTimeSeries call per variable for JLD2),
    # then use `similar` for both slots — no extra file I/O or metadata reads.
    cpu_templates, gpu_templates =
        source isa JLD2DataSource ?
            _build_templates(source, all_names, arch) :
            _build_templates(source, all_names, arch, config.grid, locations)

    gpu_frames = ntuple(_ -> Dict(k => similar(v) for (k, v) in gpu_templates), 2)
    cpu_frames = ntuple(_ -> Dict(k => similar(v) for (k, v) in cpu_templates), 2)

    reader = BufferedDataReader(source, direction,
                                Float64(config.T_start), Float64(config.T_end),
                                vel_names, var_names,
                                gpu_frames, cpu_frames,
                                1, -1, -1)

    # Prime the buffer with the first two bracketing frames at sim_t = 0
    advance_buffer!(reader, 0.0)

    return reader
end

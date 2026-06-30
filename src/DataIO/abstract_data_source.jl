"""
    AbstractDataSource

Supertype for objects that provide frame-by-frame read access to simulation output,
without requiring an intermediate copy of the data to be written to disk.

Subtypes must implement:
  - `stored_times(source)` → `Vector{Float64}` of ascending physical times in the T_start..T_end window
  - `read_frame!(cpu_field::Field, source, varname::String, idx::Int)` → fills `cpu_field` with frame `idx`
"""
abstract type AbstractDataSource end

stored_times(source::AbstractDataSource) = source.stored_times

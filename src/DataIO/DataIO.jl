module DataIO

using Oceananigans
using Oceananigans.Fields: interior, Field
using Oceananigans.OutputReaders: InMemory, FieldTimeSeries
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: on_architecture
using Oceananigans.Grids: Center, Face
using JLD2
using NCDatasets

using ..OceananigansLagrangianFilter: AbstractConfig

include("abstract_data_source.jl")
include("jld2_data_source.jl")
include("netcdf_data_source.jl")
include("buffered_data_reader.jl")

export AbstractDataSource, JLD2DataSource, NetCDFDataSource
export BufferedDataReader, create_buffered_reader
export advance_buffer!, interpolate_to_model!

end # module DataIO

# Tests for OceananigansLagrangianFilter, split across a few files by topic:
#   test_loading.jl              - the package and its exports load correctly
#   test_config_validation.jl   - OfflineFilterConfig / OnlineFilterConfig validation, esp. boundary_relaxation
#   test_initialisation.jl      - filter/map initialisation, sign flipping, and relaxation forcing
#   test_offline_integration.jl - end-to-end offline filter run vs. saved reference data
#   test_online_integration.jl  - end-to-end online filter run vs. saved reference data
using Test
using PythonCall
using Oceananigans.Units
using OceananigansLagrangianFilter

include("test_utils.jl")

include("test_loading.jl")
include("test_config_validation.jl")
include("test_initialisation.jl")
include("test_buffered_reader.jl")
include("test_offline_integration.jl")
include("test_online_integration.jl")

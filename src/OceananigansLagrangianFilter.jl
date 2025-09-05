module OceananigansLagrangianFilter

using Reexport
@reexport using Oceananigans  

# Define a supertype for all configuration objects.
abstract type AbstractConfig end

# Define submodules
include("./Utils/Utils.jl")
include("./OfflineLagrangianFilter/OfflineLagrangianFilter.jl")
include("./OnlineLagrangianFilter/OnlineLagrangianFilter.jl")

using .Utils
using .OfflineLagrangianFilter
using .OnlineLagrangianFilter

# User will have access to these functions directly
export online_test_function
export OfflineFilterConfig, run_offline_Lagrangian_filter

end # module
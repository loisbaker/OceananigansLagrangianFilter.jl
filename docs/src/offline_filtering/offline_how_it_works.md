# How it works
The offline Lagrangian filter uses many of the same helper functions as the online filter, but since the calculation doesn't need to be integrated into the user simulation, the steps are hidden and performed within the function [`run\_offline\_Lagrangian\_filter`](@ref "run_offline_Lagrangian_filter").

These steps can be seen in [`run\_offline\_lagrangian\_filter.jl`](https://github.com/loisbaker/OceananigansLagrangianFilter.jl/blob/main/src/OfflineLagrangianFilter/run_offline_lagrangian_filter.jl). 

Alternatively, the computation can be customised by directly using the helper functions rather than running [`run\_offline\_Lagrangian\_filter!`](@ref "run_offline_Lagrangian_filter!"). 

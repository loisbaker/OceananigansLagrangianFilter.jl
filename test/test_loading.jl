@testset "OceananigansLagrangianFilter loading tests" begin
# Some trivial tests to ensure that the package and key components load correctly
    # Test that key components of OceananigansLagrangianFilter are defined
    @test isdefined(OceananigansLagrangianFilter, :OfflineFilterConfig)
    @test isdefined(OceananigansLagrangianFilter, :OnlineFilterConfig)
    @test isdefined(OceananigansLagrangianFilter, :run_offline_Lagrangian_filter)
    @test isdefined(OceananigansLagrangianFilter, :create_forcing)

    # Test that Oceananigans components are also available
    @test isdefined(OceananigansLagrangianFilter, :NonhydrostaticModel)

end

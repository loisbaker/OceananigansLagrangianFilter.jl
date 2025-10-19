# Some basic tests for OceananigansLagrangianFilter
using Test
using OceananigansLagrangianFilter

@testset "OceananigansLagrangianFilter trivial tests" begin

    @test true
    # @testset "OnlineFilterConfig creation" begin
    #     using Oceananigans
    #     grid = RectilinearGrid(size=(16,16,16), x=(0,2pi), y=(0,2pi), z=(0,1))
    #     filter_config = OnlineFilterConfig(grid=grid,
    #                                        output_filename="test_output.jld2",
    #                                        var_names_to_filter=("T",),
    #                                        velocity_names=("u","v","w"),
    #                                        N=2,
    #                                        freq_c=1e-4)
    #     @test filter_config.grid === grid
    #     @test filter_config.output_filename == "test_output.jld2"
    #     @test filter_config.var_names_to_filter == ("T",)
    #     @test filter_config.velocity_names == ("u","v","w")
    #     @test filter_config.N == 2
    #     @test filter_config.freq_c == 1e-4
    # end

end
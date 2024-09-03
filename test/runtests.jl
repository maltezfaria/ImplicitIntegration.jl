using ImplicitIntegration
using Test
using Aqua

@testset "ImplicitIntegration.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(ImplicitIntegration)
    end
    # Write your tests here.
end

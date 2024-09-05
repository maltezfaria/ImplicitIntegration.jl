using Aqua
using ImplicitIntegration
using StaticArrays
using ForwardDiff
using Test

@testset "ImplicitIntegration.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(ImplicitIntegration)
    end

    @testset "Bound with IntervalArithmetic" begin
        # Utility function to test the bound function
        function test_bound(f, lc, hc, nsamples = 100^length(lc))
            N = length(lc)
            ∇f = (x) -> ForwardDiff.gradient(f, x)
            f_bnd = ImplicitIntegration.bound(f, lc, hc)
            ∇f_bnds = ImplicitIntegration.bound(∇f, lc, hc)
            for _ in 1:nsamples
                x = lc .+ (hc .- lc) .* ntuple(i -> rand(), N)
                f_bnd[1] <= f(x) <= f_bnd[2] || (return false)
                grad = ∇f(x)
                all(i -> ∇f_bnds[i][1] ≤ grad[i] ≤ ∇f_bnds[i][2], 1:N) || (return false)
            end
            return true
        end
        ## 1D
        f = x -> x[1]^3 * sin(x[1])
        lc, hc = SVector(0.0), SVector(2.0)
        @test test_bound(f, lc, hc)
        @inferred ImplicitIntegration.bound(f, lc, hc)
        ## 2D
        f = x -> x[1]^2 + x[2]^2
        lc, hc = SVector(0.0, 0.0), SVector(1.0, 1.0)
        @test test_bound(f, lc, hc)
        @inferred ImplicitIntegration.bound(f, lc, hc)
        ## 3D
        f = x -> x[1]^2 + x[2]^2 + x[3]^2
        lc, hc = SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0)
        @test test_bound(f, lc, hc)
        @inferred ImplicitIntegration.bound(f, lc, hc)
    end

    @testset "Restrict" begin
        f = (x) -> x[1]^2 + x[2]^2
        k = 1
        lc, hc = SVector(0.0, 0.0), SVector(1.0, 1.0)
        fL, fU = restrict(f, lc, hc, k)
        x̃ = SVector(0.5)
        fL(x̃)
    end
end

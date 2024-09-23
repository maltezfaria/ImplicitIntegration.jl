using Aqua
using StaticArrays
using ForwardDiff
using IntervalArithmetic
using Test
using HCubature
using ImplicitIntegration

include("aqua_tests.jl")

@testset "Bound with IntervalArithmetic" begin
    # Utility function to test the bound function
    function test_bound(f, U, nsamples = 1000)
        lc, hc = ImplicitIntegration.bounds(U)
        N = length(lc)
        ∇f = (x) -> ForwardDiff.gradient(f, x)
        f_bnd = ImplicitIntegration.bound(f, U)
        ∇f_bnds = ImplicitIntegration.bound(∇f, U)
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
    U = ImplicitIntegration.HyperRectangle(SVector(0.0), SVector(2.0))
    @test test_bound(f, U)
    @inferred ImplicitIntegration.bound(f, U)
    ## 2D
    f = x -> x[1]^2 + x[2]^2
    U = ImplicitIntegration.HyperRectangle(SVector(0.0, 0.0), SVector(1.0, 1.0))
    @test test_bound(f, U)
    @inferred ImplicitIntegration.bound(f, U)
    ## 3D
    f = x -> x[1]^2 + x[2]^2 + x[3]^2
    U = ImplicitIntegration.HyperRectangle(SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
    @test test_bound(f, U)
    @inferred ImplicitIntegration.bound(f, U)
end

@testset "Volume integrals" begin
    ## ϕ = 1
    ϕ = (x::SVector{3}) -> 1.0
    # intersection of [0,2]×[0,2]×[0,2] with {x : ϕ(x) < 0}
    U = HyperRectangle(SVector(0.0, 0.0, 0.0), SVector(2.0, 2.0, 2.0))
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 0
    # intersection of [0,2]×[0,2]×[0,2] with {x : ϕ(x) > 0}
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 2^3

    ## ϕ = z - 1
    ϕ = (x::SVector{3}) -> x[3] - 1
    # intersection of [0,2]×[0,2]×[0,2] with {x : ϕ(x) < 0}
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 2^3 / 2
    # intersection of [0,2]×[0,2]×[0,2] with {x : ϕ(x) > 0}
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 2^3 / 2

    # intersection of [0,2]×[0,2]×[0,2] with {x : ϕ(x) > 0}
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 2^3 / 2

    ## ϕ = x^2 + y^2 + z^2 - 1
    ϕ = (x::SVector{3}) -> x[1]^2 + x[2]^2 + x[3]^2 - 1
    # intersection of [0,1] × [0,1] × [0,1] with {x : ϕ(x) < 0}, i.e. a quarter
    # of a sphere
    U = HyperRectangle(SVector(0.0, 0.0, 0.0), SVector(1.1, 1.1, 1.1))
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ (4 / 3) * π / 8
end# Use quadgk as a 1D quadrature rule
quad = (f, a, b) -> hcubature(f, a, b)[1]

@testset "Segment integrals" begin
    ## ϕ = 1
    ϕ = (x::SVector{1}) -> 1.0
    # integrate on the [0,2] segment
    U = HyperRectangle(SVector(0.0), SVector(2.0))
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 2.0
    @test integrate(x -> x[1], Ω, quad) ≈ 2.0
    # empty segment
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 0
    @test integrate(x -> x[1], Ω, quad) ≈ 0

    ## ϕ = x - 1
    ϕ = (x::SVector{1}) -> x[1] - 1
    # intersection of [0,2] with {x : ϕ(x) < 0}
    U = HyperRectangle(SVector(0.0), SVector(2.0))
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 1.0
    @test integrate(x -> x[1], Ω, quad) ≈ 1 / 2
    # intersection of [0,2] with {x : ϕ(x) > 0}
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 1.0
    @test integrate(x -> x[1], Ω, quad) ≈ 2 - 1 / 2

    ## ϕ = cos(x)
    ϕ = (x::SVector{1}) -> cos(x[1])
    # intersection of [-π,π] with {x : ϕ(x) < 0}
    U = HyperRectangle(SVector(-float(π)), SVector(float(π)))
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ π
    # intersection of [-π,π] with {x : ϕ(x) < 0}
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ π

    ## ϕ = sin(x)
    ϕ = (x::SVector{1}) -> sin(x[1])
    # intersection of [-2π,2π] with {x : ϕ(x) < 0}
    U = HyperRectangle(SVector(-2 * π), SVector(2 * π))
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 2π
    # intersection of [-2π,2π] with {x : ϕ(x) < 0}
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 2π
end

@testset "Area integrals" begin
    ## ϕ = 1
    ϕ = (x::SVector{2}) -> 1.0
    # integrate on the [0,2]×[0,2] square
    U = HyperRectangle(SVector(0.0, 0.0), SVector(2.0, 2.0))
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    integrate(x -> 1.0, Ω, quad)
    @test integrate(x -> 1.0, Ω, quad) ≈ 4.0
    @test integrate(x -> x[1], Ω, quad) ≈ 4.0
    @test integrate(x -> x[2], Ω, quad) ≈ 4.0
    # empty square
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 0
    @test integrate(x -> x[1], Ω, quad) ≈ 0
    @test integrate(x -> x[2], Ω, quad) ≈ 0

    ## ϕ = y - x
    ϕ = (x::SVector{2}) -> x[2] - x[1]
    # intersection of [0,2]×[0,2] with {x : ϕ(x) < 0}
    U = HyperRectangle(SVector(0.0, 0.0), SVector(2.0, 2.0))
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 2.0
    @test integrate(x -> x[1], Ω, quad) ≈ 2^3 / 3
    @test integrate(x -> x[2], Ω, quad) ≈ 2^3 / 6

    ## ϕ = x^2 + y^2 - 1
    ϕ = (x::SVector{2}) -> x[1]^2 + x[2]^2 - 1
    # intersection of [0,1] × [0,1] with {x : ϕ(x) < 0}, i.e. quarter of a
    # circle
    U = HyperRectangle(SVector(0.0, 0.0), SVector(1.0, 1.0))
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ π / 4
    # intersection of [0,1] × [0,1] with {x : ϕ(x) > 0}
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad) ≈ 1 - π / 4
end

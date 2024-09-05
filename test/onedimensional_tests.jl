using Test
using StaticArrays: SVector
using ImplicitIntegration: HyperRectangle, ImplicitDomain, integrate
using QuadGK

# Use quadgk as a 1D quadrature rule
quad1d = (f, a, b) -> quadgk(f, a, b)[1]

@testset "Segment integrals" begin
    ## ϕ = 1
    ϕ = (x::SVector{1}) -> 1.0
    # integrate on the [0,2] segment
    U = HyperRectangle(SVector(0.0), SVector(2.0))
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad1d) ≈ 2.0
    @test integrate(x -> x[1], Ω, quad1d) ≈ 2.0
    # empty segment
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad1d) ≈ 0
    @test integrate(x -> x[1], Ω, quad1d) ≈ 0

    ## ϕ = x - 1
    ϕ = (x::SVector{1}) -> x[1] - 1
    # intersection of [0,2] with {x : ϕ(x) < 0}
    U = HyperRectangle(SVector(0.0), SVector(2.0))
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad1d) ≈ 1.0
    @test integrate(x -> x[1], Ω, quad1d) ≈ 1 / 2
    # intersection of [0,2] with {x : ϕ(x) > 0}
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad1d) ≈ 1.0
    @test integrate(x -> x[1], Ω, quad1d) ≈ 2 - 1 / 2

    ## ϕ = cos(x)
    ϕ = (x::SVector{1}) -> cos(x[1])
    # intersection of [-π,π] with {x : ϕ(x) < 0}
    U = HyperRectangle(SVector(-float(π)), SVector(float(π)))
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad1d) ≈ π
    # intersection of [-π,π] with {x : ϕ(x) < 0}
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad1d) ≈ π

    ## ϕ = sin(x)
    ϕ = (x::SVector{1}) -> sin(x[1])
    # intersection of [-2π,2π] with {x : ϕ(x) < 0}
    U = HyperRectangle(SVector(-2 * π), SVector(2 * π))
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad1d) ≈ 2π
    # intersection of [-2π,2π] with {x : ϕ(x) < 0}
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad1d) ≈ 2π
end

@testset "Area integrals" begin
    ## ϕ = 1
    ϕ = (x::SVector{2}) -> 1.0
    # integrate on the [0,2]×[0,2] square
    U = HyperRectangle(SVector(0.0, 0.0), SVector(2.0, 2.0))
    s = 1
    Ω = ImplicitDomain(U, [ϕ], [s])
    integrate(x -> 1.0, Ω, quad1d)
    @test integrate(x -> 1.0, Ω, quad1d) ≈ 4.0
    @test integrate(x -> x[1], Ω, quad1d) ≈ 4.0
    @test integrate(x -> x[2], Ω, quad1d) ≈ 4.0
    # empty square
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad1d) ≈ 0
    @test integrate(x -> x[1], Ω, quad1d) ≈ 0
    @test integrate(x -> x[2], Ω, quad1d) ≈ 0

    ## ϕ = y - x
    ϕ = (x::SVector{2}) -> x[2] - x[1]
    # intersection of [0,2]×[0,2] with {x : ϕ(x) < 0}
    U = HyperRectangle(SVector(0.0, 0.0), SVector(2.0, 2.0))
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad1d) ≈ 2.0
    @test integrate(x -> x[1], Ω, quad1d) ≈ 2^3 / 3
    @test integrate(x -> x[2], Ω, quad1d) ≈ 2^3 / 6

    ## ϕ = x^2 + y^2 - 1
    ϕ = (x::SVector{2}) -> x[1]^2 + x[2]^2 - 1
    # intersection of [0,1] × [0,1] with {x : ϕ(x) < 0}, i.e. quarter of a
    # circle
    U = HyperRectangle(SVector(0.0, 0.0), SVector(1.0, 1.0))
    s = -1
    Ω = ImplicitDomain(U, [ϕ], [s])
    @test integrate(x -> 1.0, Ω, quad1d) ≈ π / 4
end

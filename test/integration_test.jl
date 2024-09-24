using StaticArrays
using Test
using ImplicitIntegration

@testset "1D integrals" begin
    # integrate various functions on [0,2] segment
    a, b = (0.0,), (2.0,)

    ϕ = (x) -> -1.0
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 2.0
    @test integrate(x -> x[1], ϕ, a, b) ≈ 2.0

    ϕ = (x) -> 1.0
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 0
    @test integrate(x -> x[1], ϕ, a, b) ≈ 0

    ϕ = (x) -> x[1] - 1
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 1.0
    @test integrate(x -> x[1], ϕ, a, b) ≈ 1 / 2

    ϕ = (x) -> 1 - x[1]
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 1.0
    @test integrate(x -> x[1], ϕ, a, b) ≈ 2 - 1 / 2

    ϕ = (x) -> cos(x[1])
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 2 - π / 2

    ϕ = (x) -> -cos(x[1])
    @test integrate(x -> 1.0, ϕ, a, b) ≈ π / 2

    ϕ = (x) -> sin(x[1])
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 0

    ϕ = (x) -> -sin(x[1])
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 2

    ϕ = (x) -> cos(π * x[1]) # multiple roots
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 1

    @inferred integrate(x -> 1.0, ϕ, a, b)
end

@testset "2D integrals" begin
    # integrate various functions on [0,2]×[0,2] square
    a, b = (0.0, 0.0), (2.0, 2.0)

    ϕ = (x) -> -1.0
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 4.0
    @test integrate(x -> x[1], ϕ, a, b) ≈ 4.0
    @test integrate(x -> x[2], ϕ, a, b) ≈ 4.0
    @test integrate(x -> 1.0, ϕ, a, b; surface = true) ≈ 0

    ϕ = (x) -> 1.0
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 0
    @test integrate(x -> cos(x[1]), ϕ, a, b) ≈ 0
    @test integrate(x -> 1.0, ϕ, a, b; surface = true) ≈ 0

    ϕ = (x) -> x[2] - x[1]
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 2.0
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 2.0
    @test integrate(x -> x[1], ϕ, a, b) ≈ 2^3 / 3
    @test integrate(x -> x[2], ϕ, a, b) ≈ 2^3 / 6
    @test integrate(x -> 1.0, ϕ, a, b; surface = true) ≈ sqrt(8)

    ϕ = (x) -> x[1]^2 + x[2]^2 - 1
    @test integrate(x -> 1.0, ϕ, a, b) ≈ π / 4
    @test integrate(x -> 1.0, ϕ, a, b; surface = true) ≈ 2 * π / 4

    ϕ = (x) -> -(x[1]^2 + x[2]^2 - 1)
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 4 - π / 4
    @test integrate(x -> 1.0, ϕ, a, b; surface = true) ≈ 2 * π / 4

    # FIXME: type-inference fails on 1.10, but passes o 1.12. Maybe related to recursive calls in `integrate`?
    @test_broken @inferred integrate(x -> 1.0, ϕ, a, b)
end

@testset "Volume integrals" begin
    # integrate various functions on [0,2]×[0,2] square
    a, b = (0.0, 0.0, 0.0), (2.0, 2.0, 2.0)

    ϕ = (x) -> -1.0
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 2^3
    @test integrate(x -> 1.0, ϕ, a, b; surface = true) ≈ 0

    ϕ = (x) -> 1.0
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 0
    @test integrate(x -> 1.0, ϕ, a, b; surface = true) ≈ 0

    ϕ = (x) -> 1 - x[3]
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 2^3 / 2
    @test integrate(x -> 1.0, ϕ, a, b; surface = true) ≈ 2^2

    ϕ = (x) -> x[3] - 1
    @test integrate(x -> 1.0, ϕ, a, b) ≈ 2^3 / 2
    @test integrate(x -> 1.0, ϕ, a, b; surface = true) ≈ 2^2

    ϕ = (x) -> x[1]^2 + x[2]^2 + x[3]^2 - 1
    @test integrate(x -> 1.0, ϕ, a, b .+ 0.1) ≈ (4 / 3) * π / 8
    @test integrate(x -> 1.0, ϕ, a, b .+ 0.1; surface = true) ≈ 4 * π / 8

    # FIXME: type-inference fails. Maybe related to recursive calls in `integrate`?
    @test_broken @inferred integrate(x -> 1.0, ϕ, a, b .+ 0.1)
end

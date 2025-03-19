using StaticArrays
using Test
using ImplicitIntegration

order = 5
@testset "1D integrals" begin
    # integrate various functions on [0,2] segment
    a, b = (0.0,), (2.0,)

    ϕ = (x) -> -1.0
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 2.0
    @test integrate(x -> x[1], ϕ, a, b)[1] ≈ 2.0
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> x[1], Q) ≈ 2.0
    @test integrate(x -> 1.0, Q) ≈ 2.0

    ϕ = (x) -> 1.0
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 0
    @test integrate(x -> x[1], ϕ, a, b)[1] ≈ 0
    Q = quadgen(ϕ, a, b; order)[1]
    @test isempty(Q.coords) && isempty(Q.weights)

    ϕ = (x) -> x[1] - 1
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 1.0
    @test integrate(x -> x[1], ϕ, a, b)[1] ≈ 1 / 2
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> 1.0, Q) ≈ 1.0
    @test integrate(x -> x[1], Q) ≈ 1 / 2

    ϕ = (x) -> 1 - x[1]
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 1.0
    @test integrate(x -> x[1], ϕ, a, b)[1] ≈ 2 - 1 / 2
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> 1.0, Q) ≈ 1.0
    @test integrate(x -> x[1], Q) ≈ 2 - 1 / 2

    ϕ = (x) -> cos(x[1])
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 2 - π / 2
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> 1.0, Q) ≈ 2 - π / 2

    ϕ = (x) -> -cos(x[1])
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ π / 2
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> 1.0, Q) ≈ π / 2

    ϕ = (x) -> sin(x[1])
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 0
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> 1.0, Q) ≈ 0

    ϕ = (x) -> -sin(x[1])
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 2
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> 1.0, Q) ≈ 2

    ϕ = (x) -> cos(π * x[1]) # multiple roots
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 1
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> 1.0, Q) ≈ 1

    # @inferred integrate(x -> 1.0, ϕ, a, b)
    # @test_broken @inferred quadgen(ϕ, a, b; order)
end

@testset "2D integrals" begin
    # integrate various functions on [0,2]×[0,2] square
    a, b = (0.0, 0.0), (2.0, 2.0)

    ϕ = (x) -> -1.0
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 4.0
    @test integrate(x -> x[1], ϕ, a, b)[1] ≈ 4.0
    @test integrate(x -> x[2], ϕ, a, b)[1] ≈ 4.0
    @test integrate(x -> 1.0, ϕ, a, b; surface = true)[1] ≈ 0
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> 1.0, Q) ≈ 4.0
    @test integrate(x -> x[1], Q) ≈ 4.0
    @test integrate(x -> x[2], Q) ≈ 4.0
    Q = quadgen(ϕ, a, b; order, surface = true)[1]
    @test isempty(Q.coords) && isempty(Q.weights)

    ϕ = (x) -> 1.0
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 0
    @test integrate(x -> cos(x[1]), ϕ, a, b)[1] ≈ 0
    @test integrate(x -> 1.0, ϕ, a, b; surface = true)[1] ≈ 0
    Q = quadgen(ϕ, a, b; order)[1]
    @test isempty(Q.coords) && isempty(Q.weights)
    Q = quadgen(ϕ, a, b; order, surface = true)[1]
    @test isempty(Q.coords) && isempty(Q.weights)

    ϕ = (x) -> x[2] - x[1]
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 2.0
    @test integrate(x -> x[1], ϕ, a, b)[1] ≈ 2^3 / 3
    @test integrate(x -> x[2], ϕ, a, b)[1] ≈ 2^3 / 6
    @test integrate(x -> 1.0, ϕ, a, b; surface = true)[1] ≈ sqrt(8)
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> 1.0, Q)[1] ≈ 2.0
    @test integrate(x -> x[1], Q)[1] ≈ 2^3 / 3
    @test integrate(x -> x[2], Q)[1] ≈ 2^3 / 6
    Q = quadgen(ϕ, a, b; order, surface = true)[1]
    @test integrate(x -> 1.0, Q) ≈ sqrt(8)

    ϕ = (x) -> x[1]^2 + x[2]^2 - 1
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ π / 4
    @test integrate(x -> 1.0, ϕ, a, b; surface = true)[1] ≈ 2 * π / 4
    Q = quadgen(ϕ, a, b; order = 20)[1] # high-order needed for this example
    @test integrate(x -> 1.0, Q) ≈ π / 4
    Q = quadgen(ϕ, a, b; order = 20, surface = true)[1]
    @test integrate(x -> 1.0, Q) ≈ 2 * π / 4

    ϕ = (x) -> -(x[1]^2 + x[2]^2 - 1)
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 4 - π / 4
    @test integrate(x -> 1.0, ϕ, a, b; surface = true)[1] ≈ 2 * π / 4
    Q = quadgen(ϕ, a, b; order = 20)[1]
    @test integrate(x -> 1.0, Q) ≈ 4 - π / 4
    Q = quadgen(ϕ, a, b; order = 20, surface = true)[1]
    @test integrate(x -> 1.0, Q) ≈ 2 * π / 4

    # FIXME: type-inference fails on 1.10, but passes o 1.12. Maybe related to recursive calls in `integrate`?
    @test_broken @inferred integrate(x -> 1.0, ϕ, a, b)
    @test_broken @inferred quadgen(ϕ, a, b; order)
end

@testset "Volume integrals" begin
    # integrate various functions on [0,2]×[0,2] square
    a, b = (0.0, 0.0, 0.0), (2.0, 2.0, 2.0)

    ϕ = (x) -> -1.0
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 2^3
    @test integrate(x -> 1.0, ϕ, a, b; surface = true)[1] ≈ 0
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> 1.0, Q) ≈ 2^3
    Q = quadgen(ϕ, a, b; order, surface = true)[1]
    @test isempty(Q.coords) && isempty(Q.weights)

    ϕ = (x) -> 1.0
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 0
    @test integrate(x -> 1.0, ϕ, a, b; surface = true)[1] ≈ 0
    Q = quadgen(ϕ, a, b; order)[1]
    @test isempty(Q.coords) && isempty(Q.weights)

    ϕ = (x) -> 1 - x[3]
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 2^3 / 2
    @test integrate(x -> 1.0, ϕ, a, b; surface = true)[1] ≈ 2^2
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> 1.0, Q)[1] ≈ 2^3 / 2
    Q = quadgen(ϕ, a, b; order, surface = true)[1]
    @test integrate(x -> 1.0, Q) ≈ 2^2

    ϕ = (x) -> x[3] - 1
    @test integrate(x -> 1.0, ϕ, a, b)[1] ≈ 2^3 / 2
    @test integrate(x -> 1.0, ϕ, a, b; surface = true)[1] ≈ 2^2
    Q = quadgen(ϕ, a, b; order)[1]
    @test integrate(x -> 1.0, Q) ≈ 2^3 / 2
    Q = quadgen(ϕ, a, b; order, surface = true)[1]
    @test integrate(x -> 1.0, Q) ≈ 2^2

    ϕ = (x) -> x[1]^2 + x[2]^2 + x[3]^2 - 1
    @test integrate(x -> 1.0, ϕ, a, b .+ 0.1)[1] ≈ (4 / 3) * π / 8
    @test integrate(x -> 1.0, ϕ, a, b .+ 0.1; surface = true)[1] ≈ 4 * π / 8
    Q = quadgen(ϕ, a, b .+ 0.1; order = 20)[1]
    @test integrate(x -> 1.0, Q) ≈ (4 / 3) * π / 8
    Q = quadgen(ϕ, a, b .+ 0.1; order = 20, surface = true)[1]
    @test integrate(x -> 1.0, Q) ≈ 4 * π / 8

    # FIXME: type-inference fails. Maybe related to recursive calls in `integrate`?
    @test_broken @inferred integrate(x -> 1.0, ϕ, a, b .+ 0.1)
    @test_broken @inferred quadgen(ϕ, a, b .+ 0.1; order)
end

@testset "Logging" begin
    # TODO: improve testing of logger
    a, b = (0.0, 0.0, 0.0), (2.0, 2.0, 2.0)
    ϕ = (x) -> x[1]^2 + x[2]^2 + x[3]^2 - 1
    res = integrate(x -> 1.0, ϕ, a, b .+ 0.1; loginfo = true)
    @test res.val ≈ (4 / 3) * π / 8
    @test sum(res.logger.subdivisions) > 4
end

@testset "Issue 16" begin
    function dos_integer_1d_exact(E::Real, t = oneunit(E))
        x = abs(E / 2t)
        if x <= 1
            1 / sqrt(1 - x^2) / (pi * 2t)
        else
            zero(inv(oneunit(t)))
        end
    end
    ω = 0.1
    ref = dos_integer_1d_exact(ω, 0.5)
    result = integrate(
        x -> 1 / abs(2π * sinpi(2x[1])),
        x -> cospi(2x[1]) - ω,
        (0.0,),
        (1.0,);
        surface = true,
        tol = 1e-2,
    )
    @test result.val ≈ ref
end

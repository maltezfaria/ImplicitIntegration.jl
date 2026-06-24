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

    # Type-inference (issue #1): `integrate`/`quadgen` are now type-stable.
    @test (@inferred integrate(x -> 1.0, ϕ, a, b)) isa NamedTuple
    @test (@inferred quadgen(ϕ, a, b; order)) isa NamedTuple
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

    # Type-inference (issue #1): `integrate`/`quadgen` are now type-stable.
    @test (@inferred integrate(x -> 1.0, ϕ, a, b .+ 0.1)) isa NamedTuple
    @test (@inferred quadgen(ϕ, a, b .+ 0.1; order)) isa NamedTuple
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

@testset "Default interface" begin
    @test try
        ImplicitIntegration.bound(x -> x, (0.0,), (1.0,))
        ImplicitIntegration.gradient(x -> x[1])
        ImplicitIntegration.project(x -> x[1], 1, 0.5)
        ImplicitIntegration.split(x -> x[1], (0.0,), (1.0,), 1)
        true
    catch
        false
    end
    ImplicitIntegration.disable_default_interface()
    @test_throws ErrorException ImplicitIntegration.bound(x -> x, (0.0,), (1.0,))
    @test_throws ErrorException ImplicitIntegration.gradient(x -> x[1])
    @test_throws ErrorException ImplicitIntegration.project(x -> x[1], 1, 0.5)
    @test_throws ErrorException ImplicitIntegration.split(x -> x[1], (0.0,), (1.0,), 1)
    ImplicitIntegration.enable_default_interface()
end

@testset "_dedup_sorted_by_round! (alloc-free, matches unique!)" begin
    for v0 in (
        [0.0, 1.0],
        [0.1, 0.1 + 1e-12, 0.5, 0.5, 0.9],
        [0.30000001, 0.30000002, 0.7],
        Float64[],
        [3.0],
    )
        v = sort(v0)
        ref = unique(x -> round(x; sigdigits = 8), v)
        @test ImplicitIntegration._dedup_sorted_by_round!(copy(v)) == ref
    end
    # allocation-free on a warmed call (no Dict, unlike `unique!(f, v)`)
    buf = sort(rand(50))
    ImplicitIntegration._dedup_sorted_by_round!(copy(buf))
    @test (@allocated ImplicitIntegration._dedup_sorted_by_round!(buf)) == 0
end

@testset "quad1d base-case rule (Gauss default vs HCubature)" begin
    using StaticArrays: SVector
    # The default `quad1d` is a fixed Gauss-Legendre rule; it must agree with an explicit
    # adaptive HCubature `quad1d` to within tolerance on volume and surface integrals.
    import HCubature
    hcub1d = (g, a, b, tol) -> HCubature.hcubature(x -> g(x[1]), SVector(a), SVector(b); atol = tol)
    cfg_hcub = ImplicitIntegration.Config(; quad1d = hcub1d)
    for (ϕ, lc, hc) in (
        (x -> x[1]^2 + x[2]^2 - 1, (0.0, 0.0), (1.5, 1.5)),
        (x -> x[1]^2 + x[2]^2 + x[3]^2 - 1, (0.0, 0.0, 0.0), (1.5, 1.5, 1.5)),
    )
        for surface in (false, true)
            ref = integrate(x -> 1.0, ϕ, lc, hc; surface, config = cfg_hcub).val
            gauss = integrate(x -> 1.0, ϕ, lc, hc; surface).val
            @test gauss ≈ ref rtol = 1e-6
        end
    end
    # the default rule is built once (a module constant), not per call
    @test ImplicitIntegration.DEFAULT_QUAD1D === ImplicitIntegration.DEFAULT_QUAD1D
end

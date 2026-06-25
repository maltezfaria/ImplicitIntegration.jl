using StaticArrays
using Test
using ImplicitIntegration

# `integrate_threaded` splits the box into a grid of subboxes, integrates each on its own
# task and sums the results. It must agree with serial `integrate` to within the method's
# tolerance, be deterministic (independent of the number of threads), and be reentrant.
@testset "integrate_threaded" begin
    f1 = x -> 1.0

    @testset "2D agreement with serial / analytic" begin
        ϕ = x -> x[1]^2 + x[2]^2 - 1.0
        a, b = (0.0, 0.0), (1.5, 1.5)
        ser = integrate(f1, ϕ, a, b).val
        @test integrate_threaded(f1, ϕ, a, b; partition = 4).val ≈ ser rtol = 1e-6
        @test integrate_threaded(f1, ϕ, a, b; partition = 4).val ≈ π / 4 rtol = 1e-6
        # surface (perimeter)
        sser = integrate(f1, ϕ, a, b; surface = true).val
        @test integrate_threaded(f1, ϕ, a, b; partition = 4, surface = true).val ≈ sser rtol =
            1e-6
        @test integrate_threaded(f1, ϕ, a, b; partition = 4, surface = true).val ≈ 2π / 4 rtol =
            1e-6
        # a non-constant integrand (centroid numerator)
        @test integrate_threaded(x -> x[1], ϕ, a, b; partition = 5).val ≈
              integrate(x -> x[1], ϕ, a, b).val rtol = 1e-6
        # tuple partition (per-dimension subdivisions)
        @test integrate_threaded(f1, ϕ, a, b; partition = (3, 4)).val ≈ ser rtol = 1e-6
    end

    @testset "3D agreement with serial" begin
        # Off-centre sphere: no critical point (∇ϕ = 0) inside the box, so splitting never
        # triggers the low-order fallback and the threaded result matches serial tightly.
        ϕ = x -> (x[1] + 1)^2 + (x[2] + 1)^2 + (x[3] + 1)^2 - 9.0
        a, b = (0.0, 0.0, 0.0), (1.5, 1.5, 1.5)
        @test integrate_threaded(f1, ϕ, a, b; partition = 3).val ≈
              integrate(f1, ϕ, a, b).val rtol = 1e-6
        @test integrate_threaded(f1, ϕ, a, b; partition = 3, surface = true).val ≈
              integrate(f1, ϕ, a, b; surface = true).val rtol = 1e-6
    end

    @testset "determinism (independent of thread count)" begin
        ϕ = x -> x[1]^2 + x[2]^2 - 1.0
        a, b = (0.0, 0.0), (1.5, 1.5)
        # results are summed in fixed grid order, so repeated calls are bit-for-bit equal
        r1 = integrate_threaded(f1, ϕ, a, b; partition = 7).val
        r2 = integrate_threaded(f1, ϕ, a, b; partition = 7).val
        @test r1 === r2
    end

    @testset "type stability and API" begin
        ϕ = x -> x[1]^2 + x[2]^2 - 1.0
        a, b = (0.0, 0.0), (1.5, 1.5)
        @test integrate_threaded(f1, ϕ, a, b; partition = 2).val isa Float64
        # loginfo is unsupported and ignored (with a warning), logger is nothing
        res = @test_logs (:warn,) integrate_threaded(f1, ϕ, a, b; partition = 2, loginfo = true)
        @test res.logger === nothing
        @test_throws ArgumentError integrate_threaded(f1, ϕ, a, b; partition = 0)
    end
end

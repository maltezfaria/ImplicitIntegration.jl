using ImplicitIntegration
using Test
using StaticArrays
using DynamicPolynomials

@testset "Evaluation" begin
    @testset "1D" begin
        n = 2 # degree of the Bernstein polynomial
        lc = SVector(0.0)
        hc = SVector(1.0)
        p = ImplicitIntegration.BernsteinPolynomial([1.0, 0, 0], lc, hc) # p(x) = (1 - x)^2
        x = SVector(0.1)
        @test p(x) ≈ (1 - x[1])^2
        p = ImplicitIntegration.BernsteinPolynomial([0.0, 1.0, 0], lc, hc) # p(x) = 2x(1-x)
        @test p(x) ≈ 2 * x[1] * (1 - x[1])
    end
    @testset "2D" begin
        n = 2 # degree of the Bernstein polynomial
        lc = SVector(0.0, 0.0)
        hc = SVector(1.0, 1.0)
        c = zeros(n + 1, n + 1)
        c[1, 1] = 1.0 # (1-x₁)^2*(1-x₂)^2
        p = ImplicitIntegration.BernsteinPolynomial(c, lc, hc)
        x = SVector(0.1, 0.2)
        @test p(x) ≈ (1 - x[1])^2 * (1 - x[2])^2
        fill!(c, 0.0)
        c[2, 1] = 1.0 # 2x₁(1-x₁)(1-x₂)^2
        @test p(x) ≈ 2 * x[1] * (1 - x[1]) * (1 - x[2])^2
        c[1, 2] = 2.0 # += 4x₂(1-x₂)(1-x₁)^2
        @test p(x) ≈
              2 * x[1] * (1 - x[1]) * (1 - x[2])^2 + 4 * x[2] * (1 - x[2]) * (1 - x[1])^2
    end
end

@testset "Derivative" begin
    @testset "1D" begin
        n = 2 # degree of the Bernstein polynomial
        lc = SVector(0.0)
        hc = SVector(1.0)
        p = ImplicitIntegration.BernsteinPolynomial([1.0, 0, 0], lc, hc) # p(x) = (1 - x)^2
        dp = ImplicitIntegration.derivative(p, 1)
        @test dp(SVector(0.1)) ≈ -2 * (1 - 0.1)
    end
    @testset "2D" begin
        n = 2 # degree of the Bernstein polynomial
        lc = SVector(0.0, 0.0)
        hc = SVector(1.0, 1.0)
        c = zeros(n + 1, n + 1)
        c[1, 1] = 1.0 # (1-x₁)^2*(1-x₂)^2
        p = ImplicitIntegration.BernsteinPolynomial(c, lc, hc)
        dp = ImplicitIntegration.derivative(p, 1)
        @test dp(SVector(0.1, 0.2)) ≈ -2 * (1 - 0.1) * (1 - 0.2)^2
        dp = ImplicitIntegration.derivative(p, 2)
        @test dp(SVector(0.1, 0.2)) ≈ -2 * (1 - 0.1)^2 * (1 - 0.2)
        ∇p = ImplicitIntegration.gradient(p)
        x = SVector(0.1, 0.2)
        @test ∇p(x) ≈
              SVector(-2 * (1 - x[1]) * (1 - x[2])^2, -2 * (1 - x[1])^2 * (1 - x[2]))
    end
end

@testset "Split" begin
    @testset "1D" begin
        lc = SVector(1.1)
        hc = SVector(2.3)
        p = ImplicitIntegration.BernsteinPolynomial(rand(3), lc, hc)
        α = 0.6
        p1, p2 = ImplicitIntegration.split(p, 1, α)
        @test ImplicitIntegration.low_corner(p1) == lc
        @test ImplicitIntegration.high_corner(p1) == lc + α * (hc - lc)
        @test ImplicitIntegration.low_corner(p2) == lc + α * (hc - lc)
        @test ImplicitIntegration.high_corner(p2) == hc
        x1 = lc + α / 2 * (hc - lc)
        @test p1(x1) ≈ p(x1)
        x2 = lc + (1 - α) / 2 * (hc - lc)
        @test p2(x2) ≈ p(x2)
    end
    @testset "2D" begin
        lc = SVector(1.1, 2.2)
        hc = SVector(2.3, 3.5)
        c = rand(3, 3)
        p = ImplicitIntegration.BernsteinPolynomial(c, lc, hc)
        d = 1
        α = 0.6
        p1, p2 = ImplicitIntegration.split(p, d, α)
        @test ImplicitIntegration.low_corner(p1) == lc
        @test ImplicitIntegration.high_corner(p1)[d] == lc[d] + α * (hc[d] - lc[d])
        @test ImplicitIntegration.low_corner(p2)[d] == lc[d] + α * (hc[d] - lc[d])
        @test ImplicitIntegration.high_corner(p2) == hc
        x1 = lc + α / 2 * (hc - lc)
        @test p1(x1) ≈ p(x1)
        x2 = lc + (1 - α) / 2 * (hc - lc)
        @test p2(x2) ≈ p(x2)
    end
end

@testset "Restrict" begin
    @testset "2D" begin
        n = 2 # degree of the Bernstein polynomial
        lc = SVector(0.0, 0.0)
        hc = SVector(1.0, 1.0)
        c = zeros(n + 1, n + 1)
        c[1, 1] = 1.0 # (1-x₁)^2*(1-x₂)^2
        p = ImplicitIntegration.BernsteinPolynomial(c, lc, hc)
        d = 1
        p_low = ImplicitIntegration.lower_restrict(p, d) # restrict to x₁ = 0
        @test all(p_low(x) ≈ (1 - x[1])^2 for x in rand(SVector{1}, 10))
        p_high = ImplicitIntegration.upper_restrict(p, d) # restrict to x₁ = 1
        @test all(p_high(x) ≈ 0 for x in rand(SVector{1}, 10))
        #
        d = 2
        p_low = ImplicitIntegration.lower_restrict(p, d) # restrict to x₂ = 0
        @test all(p_low(x) ≈ (1 - x[1])^2 for x in rand(SVector{1}, 10))
        p_high = ImplicitIntegration.upper_restrict(p, d) # restrict to x₂ = 1
        @test all(p_high(x) ≈ 0 for x in rand(SVector{1}, 10))
    end
end

@testset "Polynomial" begin
    coefs = (Dict((2, 0) => 1.0, (0, 2) => 1.0, (0, 0) => -1.0))
    p = ImplicitIntegration.Polynomial(coefs) # p(x) = x₁² + x₂² - 1
    ϕ = (x) -> x[1]^2 + x[2]^2 - 1
    lb, ub = (-2.0, -2.0), (2.0, 2.0)
    b = ImplicitIntegration.BernsteinPolynomial(p, lb, ub)
    @test all(b(x) ≈ ϕ(x) for x in rand(SVector{2}, 1000))
end

@testset "Integrate" begin
    ImplicitIntegration.disable_default_interface()
    coefs = (Dict((2, 0) => 1.0, (0, 2) => 1.0, (0, 0) => -1.0))
    p = ImplicitIntegration.Polynomial(coefs) # p(x) = x₁² + x₂² - 1
    ϕ = (x) -> x[1]^2 + x[2]^2 - 1
    lb, ub = (-2.0, -2.0), (2.0, 2.0)
    b = ImplicitIntegration.BernsteinPolynomial(p, lb, ub)
    @test integrate(x -> 1.0, b, lb, ub)[1] ≈ π
    Q = quadgen(b, lb, ub; order = 20)[1]
    @test integrate(x -> 1.0, Q) ≈ π
    ImplicitIntegration.enable_default_interface()
end

@testset "DynamicPolynomials" begin
    ImplicitIntegration.disable_default_interface()
    lb, ub = (-2.0, -2.0), (2.0, 2.0)
    @polyvar x y
    p = x^2 + y^2 - 1.0
    b = ImplicitIntegration.BernsteinPolynomial(p, (-2.0, -2.0), (2.0, 2.0))
    @test all(b(x) ≈ p(x) for x in rand(SVector{2}, 1000))
    @test integrate(x -> 1.0, b, lb, ub).val ≈ π
    ImplicitIntegration.enable_default_interface()
end

# @testset "Interpolation" begin
#     @testset "1D" begin
#         f = (x) -> (1 - x[1])^2 + x[1]^4
#         degree = (100,)
#         lb = SVector(0.0)
#         ub = SVector(1.0)
#         x = ImplicitIntegration.uniform_points(degree, lb, ub)
#         p = ImplicitIntegration.BernsteinPolynomial(f.(x), x, lb, ub)
#         @test all(f(x) ≈ p(x) for x in rand(SVector{1}, 1000))
#         V = ImplicitIntegration.vandermonde_matrix(degree, x)
#         vals = f.(x)
#         p = ImplicitIntegration.BernsteinPolynomial(V \ vals, SVector(0.0), SVector(1.0))
#         @test all(f(x) ≈ p(x) for x in rand(SVector{1}, 1000))
#     end
#     @testset "2D" begin
#         f = (x) -> (1 - x[1])^2 + x[1]^4 + x[2]^5 * x[1]^3
#         degree = (4, 5)
#         lc = SVector(0.0, 0.0)
#         hc = SVector(1.0, 1.0)
#         x = ImplicitIntegration.uniform_points(degree, lc, hc)
#         V = ImplicitIntegration.vandermonde_matrix(degree, x)
#         vals = f.(x)
#         c = reshape(V \ vec(vals), degree .+ 1)
#         p = ImplicitIntegration.BernsteinPolynomial(c, lc, hc)
#         @test all(f(x) ≈ p(x) for x in rand(SVector{2}, 1000))
#     end
# end

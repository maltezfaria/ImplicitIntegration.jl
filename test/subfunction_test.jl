using ImplicitIntegration
using StaticArrays
using IntervalArithmetic
using ForwardDiff
using Test

f = (x::SVector{3}) -> x[1]^2 + x[2]^2 + x[3]^2
idxs = @SVector [1, 3]
vals = @SVector rand(2)
f̃ = ImplicitIntegration.SubFunction{1}(f, idxs, vals)
g = (x) -> f(SVector(vals[1], x[1], vals[2]))
x̃ = @SVector rand(1)
x = ImplicitIntegration.multi_insert(x̃, idxs, vals)
@test x == SVector(vals[1], x̃[1], vals[2])
@test @inferred(f̃(x̃)) == f(x) == g(x̃)

# bound
U = ImplicitIntegration.HyperRectangle(SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
Ũ = ImplicitIntegration.HyperRectangle(SVector(0.0), SVector(1.0))
@test ImplicitIntegration.bound(f̃, Ũ) isa Tuple{Float64,Float64}
@test @inferred(ImplicitIntegration.bound(f̃, Ũ)) == ImplicitIntegration.bound(g, Ũ)

# gradient
∇g = ForwardDiff.gradient(g, x̃)
∇f̃ = @inferred(ImplicitIntegration.gradient(f̃, x̃))
@test ∇f̃ == ∇g

# bound_gradient
intrv = interval(0.0, 1.0)
∇g_bnd = IntervalArithmetic.bounds.(ForwardDiff.gradient(g, SVector(intrv)))
∇f̃_bnd = @inferred(ImplicitIntegration.bound_gradient(f̃, Ũ))
@test ∇f̃_bnd == ∇g_bnd

# try a different restriction
idxs = @SVector [2]
vals = @SVector rand(1)
f̃ = ImplicitIntegration.SubFunction{2}(f, idxs, vals)
g = (x) -> f(SVector(x[1], vals[1], x[2]))
x̃ = @SVector rand(2)
x = ImplicitIntegration.multi_insert(x̃, idxs, vals)
@test x == SVector(x̃[1], vals[1], x̃[2])
@test @inferred(f̃(x̃)) == f(x) == g(x̃)

# bound
U = ImplicitIntegration.HyperRectangle(SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
Ũ = ImplicitIntegration.remove_dimension(U, 2)
@test ImplicitIntegration.bound(f̃, Ũ) isa Tuple{Float64,Float64}
@test @inferred(ImplicitIntegration.bound(f̃, Ũ)) == ImplicitIntegration.bound(g, Ũ)

# gradient
∇g = ForwardDiff.gradient(g, x̃)
∇f̃ = @inferred(ImplicitIntegration.gradient(f̃, x̃))
@test ∇f̃ == ∇g

# bound_gradient
intrv = (interval(0.0, 1.0), interval(0.0, 1.0))
∇g_bnd = IntervalArithmetic.bounds.(ForwardDiff.gradient(g, SVector(intrv)))
∇f̃_bnd = @inferred(ImplicitIntegration.bound_gradient(f̃, Ũ))
@test ∇f̃_bnd == ∇g_bnd

# try the "trivial" case
idxs = SVector{0,Int}()
vals = SVector{0,Float64}()
f̃ = ImplicitIntegration.SubFunction{3}(f, idxs, vals)
x = @SVector rand(3)
@test @inferred(f̃(x)) == f(x)

# bound
U = ImplicitIntegration.HyperRectangle(SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
@test ImplicitIntegration.bound(f̃, U) isa Tuple{Float64,Float64}
@test @inferred(ImplicitIntegration.bound(f̃, U)) == ImplicitIntegration.bound(f, U)

# gradient
∇f̃ = @inferred(ImplicitIntegration.gradient(f̃, x))
@test ∇f̃ == ForwardDiff.gradient(f, x)

# bound_gradient
intrv = (interval(0.0, 1.0), interval(0.0, 1.0), interval(0.0, 1.0))
∇g_bnd = IntervalArithmetic.bounds.(ForwardDiff.gradient(f, SVector(intrv)))
∇f̃_bnd = @inferred(ImplicitIntegration.bound_gradient(f̃, U))
@test ∇f̃_bnd == ∇g_bnd

#=
The following methods constitute the interace required by any function that is
passed to the `integrate` function. By default we use `ForwardDiff` to compute
gradients and `IntervalArithmetic` to compute bounds. However, the user can
overload these to use a custom implementation.
=#

"""
    bound(f, lc, hc) --> (lb, ub)

Return a lower and upper bound for the function `f : U → ℝ` valid for all `x ∈ U`.
"""
function bound(f, lc, hc)
    N = length(lc)
    I = ntuple(i -> IntervalArithmetic.interval(lc[i], hc[i]), N) |> SVector
    return IntervalArithmetic.bounds.(f(I))
end
bound(f, rec) = bound(f, bounds(rec)...)

"""
    gradient(f, x)

Compute the gradient of a function `f : ℝᵈ → ℝ` at point `x ∈ ℝᵈ`.
"""
function gradient(f, x)
    return ForwardDiff.gradient(f, x)
end

"""
    bound_gradient(f, lc, hc) --> bnds

Compute a lower and upper bound for the gradient of a function `f : U → ℝ` valid
for all `x ∈ U` in the sense that `bnds[i][1] ≤ ∂f/∂xᵢ(x) ≤ bnds[i][2]`.
"""
function bound_gradient(f, lc, hc)
    N = length(lc)
    ∇f = x -> ForwardDiff.gradient(f, x)
    I = ntuple(i -> IntervalArithmetic.interval(lc[i], hc[i]), N) |> SVector
    return IntervalArithmetic.bounds.(∇f(I))
end
bound_gradient(f, rec::HyperRectangle) = bound_gradient(f, bounds(rec)...)

"""
    restrict(f, lc, hc, k) --> (fl, fu)

Given a function `f : U → ℝ`, where `U = {lc ≤ x ≤ hc}`, return two functions `fl` and `fu`
representing the restriction of `f` to upper and lower faces of `U` along the `k`-th
dimension.

Note that the returned functions are required to also implement the interface defined here.
"""
function restrict(f, lc, hc, k)
    M = length(lc)
    val_lower = lc[k]
    val_upper = hc[k]
    fl = SubFunction{M - 1}(f, SVector(k), SVector(val_lower))
    fu = SubFunction{M - 1}(f, SVector(k), SVector(val_upper))
    return fl, fu
end
restrict(f, rec, k) = restrict(f, bounds(rec)..., k)

## Heuristic "bounds" based on the function values at the corners of the rectangle

function heuristic_bound(f, lc, hc, n = 10)
    iters = ntuple(i -> range(lc[i], hc[i], n), length(lc))
    extrema(Iterators.product(iters...)) do x
        return f(SVector(x))
    end
end

function heuristic_bound_gradient(f, lc, hc, n = 10)
    N = length(lc)
    lbnds = svector(i -> Inf, N)
    ubnds = svector(i -> -Inf, N)
    iters = ntuple(i -> range(lc[i], hc[i], n), N)
    for x in Iterators.product(iters...)
        v = gradient(f, SVector(x))
        lbnds = min.(lbnds, v)
        ubnds = max.(ubnds, v)
    end
    return svector(i -> SVector(lbnds[i], ubnds[i]), N)
end

function use_heuristic_bounds(F::Type, n = 10)
    @eval begin
        function ImplicitIntegration.bound(f::$F, lc, hc)
            return ImplicitIntegration.heuristic_bound(f, lc, hc, $n)
        end

        function ImplicitIntegration.bound_gradient(f::$F, lc, hc)
            return ImplicitIntegration.heuristic_bound_gradient(f, lc, hc, $n)
        end
    end
end
use_heuristic_bounds(f, args...) = use_heuristic_bounds(typeof(f), args...)

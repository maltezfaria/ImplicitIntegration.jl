#=
The following methods constitute the interace required by any function that is passed to the
`integrate` function. By default we use `ForwardDiff` to compute gradients and
`IntervalArithmetic` to compute bounds. However, the user can overload these to use a custom
implementation.
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
     restrict(f, k, v)

Given a function `f : ℝᵈ → ℝ`, a value `v ∈ ℝ` and an integer `1 ≤ k ≤ d`, return the
function `f̃ : ℝᵈ⁻¹ → ℝ` defined by restricting `f` to the value `v` along dimension `d`;
i.e. `f̃(x) = f(x₁, ..., x_{k-1}, v, x_{k}, ..., x_d)`.

!!! note

    The returned type should also implement the interface methods `gradient`, `bound`,
    and `bound_gradient`.
"""
restrict(f, k, v) = (x) -> f(insert(x, k, v))

## Heuristic "bounds" based on the function values at the corners of the rectangle

function heuristic_bound(f, lc, hc, n = 10)
    iters = ntuple(i -> range(lc[i], hc[i], n), length(lc))
    extrema(Iterators.product(iters...)) do x
        return f(SVector(x))
    end
end

function heuristic_bound_gradient(f, lc, hc, n = 10)
    N = length(lc)
    lbnds = ntuple(i -> Inf, N) |> SVector
    ubnds = ntuple(i -> -Inf, N) |> SVector
    iters = ntuple(i -> range(lc[i], hc[i], n), N)
    for x in Iterators.product(iters...)
        v = gradient(f, SVector(x))
        lbnds = min.(lbnds, v)
        ubnds = max.(ubnds, v)
    end
    bnds = ntuple(i -> SVector(lbnds[i], ubnds[i]), N)
    return SVector(bnds)
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

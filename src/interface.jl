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
    gradient(f)

Compute the gradient function `f : ℝᵈ → ℝ`. The returned function takes a vector `x ∈ ℝᵈ`
and returns the gradient `∇f(x) ∈ ℝᵈ`.
"""
function gradient(f)
    return x -> ForwardDiff.gradient(f, x)
end

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

"""
    use_heuristic_bounds(F::Type, n = 10)

Overload interface methods for type `F` to use a sample-based heuristic when bounding
functions of type `F` and its gradient.

The bounds are obtained by sampling the function (or its gradient) on an `n × ... × n` grid
and taking the maximum/minimum of the attained values. This is obviously a heuristic, and
may fail in practice.
"""
function use_heuristic_bounds(F::Type, n = 10)
    @eval begin
        function ImplicitIntegration.bound(f::$F, lc, hc)
            return ImplicitIntegration.heuristic_bound(f, lc, hc, $n)
        end
    end
end
use_heuristic_bounds(f, args...) = use_heuristic_bounds(typeof(f), args...)

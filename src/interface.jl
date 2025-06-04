#=
The following methods constitute the interace required by any function that is passed to the
`integrate` function. By default we use `ForwardDiff` to compute gradients and
`IntervalArithmetic` to compute bounds. However, the user can overload these to use a custom
implementation.
=#

# easily disable the default interface for e.g. testing that all needed methods are
# implemented
const ALLOW_DEFAULT_INTERFACE = Ref(true)

function disable_default_interface()
    return ALLOW_DEFAULT_INTERFACE[] = false
end

function enable_default_interface()
    return ALLOW_DEFAULT_INTERFACE[] = true
end

"""
    bound(f, lc, hc) --> (lb, ub)

Return a lower and upper bound for the function `f : U → ℝ` valid for all `x ∈ U`.

    bound(∇f, lc, hc) --> bnds

Compute a lower and upper bound for the gradient of a function `f : U → ℝ` valid
for all `x ∈ U` in the sense that `bnds[i][1] ≤ ∂f/∂xᵢ(x) ≤ bnds[i][2]`.
"""
function bound(f, lc, hc)
    ALLOW_DEFAULT_INTERFACE[] || error(
        "Default interface is disabled, please implement the `bound` method for your function type.",
    )
    N = length(lc)
    I = ntuple(i -> IntervalArithmetic.interval(lc[i], hc[i]), N) |> SVector
    return IntervalArithmetic.bounds.(f(I))
end
bound(f, rec::HyperRectangle) = bound(f, bounds(rec)...)

"""
    gradient(f)

Compute the gradient function `f : ℝᵈ → ℝ`. The returned function takes a vector `x ∈ ℝᵈ`
and returns the gradient `∇f(x) ∈ ℝᵈ`.
"""
function gradient(f)
    ALLOW_DEFAULT_INTERFACE[] || error(
        "Default interface is disabled, please implement the `gradient` method for your function type.",
    )
    return x -> ForwardDiff.gradient(f, x)
end

"""
     project(f, k, v)

Given a function `f : ℝᵈ → ℝ`, a value `v ∈ ℝ` and an integer `1 ≤ k ≤ d`, return the
function `f̃ : ℝᵈ⁻¹ → ℝ` defined by projecting `f` onto the hyperplane `xₖ = v`; i.e. `f̃(x) = f(x₁, ..., x_{k-1}, v, x_{k}, ..., x_d)`.

!!! note

    The returned type should also implement the interface methods `gradient`, `bound`.
"""
function project(f, k, v)
    ALLOW_DEFAULT_INTERFACE[] || error(
        "Default interface is disabled, please implement the `project` method for your function type.",
    )
    return (x) -> f(insert(x, k, v))
end

"""
    split(f, lb, ub)

Given a function `f : ℝᵈ → ℝ` and lower and upper bounds `lb, ub ∈ ℝᵈ`, return the
restriction of `f` to the hyperrectangle defined by `lb` and `ub`.

By default this function simply returns `f`, but it computing sharper bounds on the
restricted function requires a more sophisticated implementation.
"""
function split(f, lb, ub, dir)
    ALLOW_DEFAULT_INTERFACE[] || error(
        "Default interface is disabled, please implement the `split` method for your function type.",
    )
    return f, f
end
function split(f, rec::HyperRectangle, dir)
    lc, hc = bounds(rec)
    return split(f, lc, hc, dir)
end

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

#=
The following methods constitute the interace required by any function that is
passed to the `integrate` function. By default we use `ForwardDiff` to compute
gradients and `IntervalArithmetic` to compute bounds. However, the user can
overload these to use a custom implementation.
=#

"""
    bound(f, U::HyperRectangle) --> (lb, ub)

Return a lower and upper bound for the function `f : U → ℝ` valid for all `x ∈
U`.
"""
function bound(f, rec)
    lc, hc = bounds(rec)
    I = ntuple(i -> IntervalArithmetic.interval(lc[i], hc[i]), ndims(rec)) |> SVector
    return IntervalArithmetic.bounds.(f(I))
end

"""
    gradient(f, x)

Compute the gradient of a function `f : ℝᵈ → ℝ` at point `x ∈ ℝᵈ`.
"""
function gradient(f, x)
    return ForwardDiff.gradient(f, x)
end

"""
    bound_gradient(f, U::HyperRectangle) --> (lbs, ubs)

Compute a lower and upper bound for the gradient of a function `f : U → ℝ` valid
for all `x ∈ U` in the sense that `lbs[i] ≤ ∂f/∂xᵢ(x) ≤ ubs[i]`.
"""
function bound_gradient(f, rec)
    ∇f = x -> gradient(f, x)
    return bound(∇f, rec)
end

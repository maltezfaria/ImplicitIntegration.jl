#=
The following methods constitute the interace required by any function that is
passed to the `integrate` function. By default we use `ForwardDiff` to compute
gradients and `IntervalArithmetic` to compute bounds. However, the user can
overload these to use a custom implementation.
=#

"""
    bound(f, lc, hc) --> (lb, ub)

Return a lower and upper bound for the function `f : U → ℝ` valid for all `x ∈
U`.
"""
function bound(f, lc, hc)
    I = ntuple(i -> IntervalArithmetic.interval(lc[i], hc[i]), ndims(rec)) |> SVector
    return IntervalArithmetic.bounds.(f(I))
end
bound(f, rec) = bound(f, bounds(rec))

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
    ∇f = x -> gradient(f, x)
    return bound(∇f, lc, hc)
end
bound_gradient(f, rec) = bound_gradient(f, bounds(rec))

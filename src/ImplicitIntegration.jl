module ImplicitIntegration

import IntervalArithmetic
import ForwardDiff
import Roots
import HCubature
import StaticArrays: SVector, MVector, insert, deleteat, setindex, popfirst

# FIXME: understand the patch below and move it upstream if necessary
IntervalArithmetic.Interval{T}(x::T) where {T<:Real} = IntervalArithmetic.interval(x)

include("utils.jl")
include("hyperrectangle.jl")
include("subfunction.jl")
include("integration.jl")

"""
    bound(f, U::HyperRectangle) --> (lb, ub)

Return a lower and upper bound for the function `f : U → ℝ` valid for all `x ∈
U`.
"""
function bound(f, rec::HyperRectangle{N}) where {N}
    lc, hc = bounds(rec)
    I = ntuple(i -> IntervalArithmetic.interval(lc[i], hc[i]), N) |> SVector
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
function bound_gradient(f, rec::HyperRectangle{N}) where {N}
    ∇f = x -> gradient(f, x)
    return bound(∇f, rec)
end

export integrate

end # module

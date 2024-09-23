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
include("interface.jl")
include("subfunction.jl")
include("integration.jl")

export integrate

end # module

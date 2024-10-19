module ImplicitIntegration

using LinearAlgebra

import ForwardDiff
import IntervalArithmetic
import Roots
import HCubature
import StaticArrays: SVector, MVector, insert, deleteat, setindex, popfirst, pop, push

# FIXME: understand the patch below and move it upstream if necessary
IntervalArithmetic.Interval{T}(x) where {T<:Real} = IntervalArithmetic.interval(x)

include("utils.jl")
include("hyperrectangle.jl")
include("interface.jl")
include("subfunction.jl")
include("integration.jl")
include("quadgen.jl")

export integrate, quadgen

end # module

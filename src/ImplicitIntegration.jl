module ImplicitIntegration

using LinearAlgebra

import ForwardDiff
import IntervalArithmetic
import Roots
import HCubature
import StaticArrays: SVector, MVector, insert, deleteat, setindex, popfirst, pop, push
import Memoize

include("hyperrectangle.jl")
include("treenode.jl")
include("interface.jl")
include("integration.jl")
include("quadgen.jl")
include("threaded.jl")
include("bernstein.jl")

export integrate, integrate_threaded, quadgen

end # module

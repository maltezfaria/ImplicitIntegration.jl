module ImplicitIntegration

import IntervalArithmetic
import ForwardDiff
import Roots
import HCubature
import StaticArrays: SVector, insert, deleteat, setindex

# FIXME: understand the patch below and move it upstream if necessary
IntervalArithmetic.Interval{T}(x::T) where {T} = IntervalArithmetic.interval(x)

include("hyperrectangle.jl")

"""
    bound(f, rec::HyperRectangle)

Return a lower and upper bound for the function `f` on the hypercube `rec`.

When `f` is a scalar function, return a lower and upper bound. When `f` is
vector valued, return a vector of lower and upper bounds.
"""
function bound(f, rec::HyperRectangle{N}) where {N}
    lc, hc = bounds(rec)
    I = ntuple(i -> IntervalArithmetic.interval(lc[i], hc[i]), N) |> SVector
    return IntervalArithmetic.bounds.(f(I))
end

function restrict(f, U::HyperRectangle{N}, k) where {N}
    @assert 1 ≤ k ≤ N
    lc, hc = bounds(U)
    fl = x̃ -> f(insert(x̃, k, lc[k]))
    fu = x̃ -> f(insert(x̃, k, hc[k]))
    return fl, fu
end

"""
    struct ImplicitDomain

Structure used to represent a domain `V` as the intersection of several
implicitly defined domains `Vᵢ`. The domains `Vᵢ` are defined through a set of
`functions` and `signs` as follows:

- `Vᵢ = {x ∈ U : ϕᵢ(x) > 0}` if `sᵢ = 1`.
- `Vᵢ = {x ∈ U : ϕᵢ(x) < 0}` if `sᵢ = -1`.
- `Vᵢ = {x ∈ U : ϕᵢ(x) ≠ 0}` if `sᵢ = 0`.
"""
struct ImplicitDomain{N,T}
    box::HyperRectangle{N,T}
    functions::Vector{<:Function}
    flags::Vector{Int}
end

## Algorithm 1 of Saye 2015
"""
    integrand_eval(f, V::ImplicitDomain{N}, dir, quad1d)

Return a function `f̃ : ℝᴺ⁻¹ -> ℝ` that approximates the one-dimensional
integral over `I(x̃) ⊂ ℝ` of the function `t -> f(insert(x̃, k, t))`, where the
integration domain `I` is defined as `I(x̃) = {t ∈ [a,b] :
sᵢ*ϕᵢ(insert(̃x,k,t) ≥ 0 ∀ (ϕᵢ,sᵢ) ∈ V)}`.
"""
function integrand_eval(f, V::ImplicitDomain{N}, dir::Int, quad) where {N}
    xl, xu = bounds(V.box)
    a, b = xl[dir], xu[dir]
    f̃ = (x̃) -> begin
        @assert length(x̃) == N - 1 "dimension of x̃ must be $(N - 1)"
        # compute the connected components
        bnds = [a, b]
        for ϕᵢ in V.functions
            g = (t) -> ϕᵢ(insert(x̃, dir, t))
            if N == 1
                # we need to find all roots of g in [a,b]
                append!(bnds, Roots.find_zeros(g, a, b))
            else
                # we know that g is monotonic since it corresponds to a
                # height-direction, so at most a single root exists.
                # TODO: look up and test Brent's method for this case
                g(a) * g(b) > 0 && continue
                push!(bnds, Roots.find_zero(g, (a, b)))
            end
        end
        sort!(bnds)
        I = 0.0 # FIXME: infer the appropriate type for this accumulator
        # compute the integral
        for i in 1:(length(bnds)-1) # loop over each segment
            rᵢ, rᵢ₊₁ = bnds[i], bnds[i+1]
            L = rᵢ₊₁ - rᵢ
            # decide if the segment is inside the domain
            xc = insert(x̃, dir, (rᵢ + rᵢ₊₁) / 2)
            skip = false
            for (ϕᵢ, sᵢ) in zip(V.functions, V.flags)
                sᵢ == 0 && continue # avoid evaluation of ϕᵢ(xc) when possible
                if sᵢ * ϕᵢ(xc) < 0
                    skip = true
                    break
                end
            end
            skip && continue
            # add the contribution of the segment by performing a 1D quadrature.
            I += quad(SVector(rᵢ), SVector(rᵢ₊₁)) do (t,)
                x = insert(x̃, dir, t)
                return f(x)
            end
        end
        return I
    end
    return f̃
end

function integrate(f, Ω::ImplicitDomain{DIM,T}, quad) where {DIM,T}
    if DIM == 1
        # Base case of the dimensional recursion
        f̃ = integrand_eval(f, Ω, 1, quad)
        x̃ = SVector{0,T}() # TODO: this is a hack to make `integrand_eval` work for N=1
        return f̃(x̃)
    else
        # First we prune the domain by removing functions that are not relevant
        U = Ω.box
        xl, xu = bounds(U)
        xc = (xl + xu) / 2
        ϕ, s = Ω.functions, Ω.flags # before prunning
        n = length(ϕ)
        idxs_keep = Int[]
        for i in 1:n
            lb, ub = bound(ϕ[i], U)
            lb > 0 && s[i] < 0 && (return 0.0) # {x ∈ U : ϕᵢ < 0} is empty since ϕᵢ > 0
            ub < 0 && s[i] > 0 && (return 0.0) # {x ∈ U : ϕᵢ > 0} is empty since ϕᵢ < 0
            lb * ub ≤ 0 && push!(idxs_keep, i) # sign change in U, so keep ϕᵢ
            # NOTE: things like lb > 0 and s[i] > 0 can be pruned since U ∩ {x :
            # ϕᵢ > 0} = U
        end
        # Pruned version
        ϕ, s = ϕ[idxs_keep], s[idxs_keep]
        n    = length(idxs_keep)
        # tensor product quadrature when prunning removes all ϕᵢ
        if n == 0
            return quad(f, xl, xu)
        end
        # Domain is neither empty nor full, so we need to either split it or
        # find a height direction...
        ∇ϕ₁ = (x) -> ForwardDiff.gradient(ϕ[1], x)
        k = argmax(∇ϕ₁(xc))
        ϕ̃, s̃ = Function[], empty(s)
        for i in 1:n
            ∇ϕᵢ = (x) -> ForwardDiff.gradient(ϕ[i], x)
            ∇ϕᵢ_bnds = bound(∇ϕᵢ, U)
            den = sum(∇ϕᵢ_bnds) do (lb, ub)
                return max(abs(lb), abs(ub))^2
            end |> sqrt # max over U of |∇ϕᵢ|
            lb, ub = ∇ϕᵢ_bnds[k]
            min_qual = 0
            qual = den == 0 ? 1.0 : lb * ub > 0 ? min(abs(lb), abs(ub)) / den : 0.0 # quality factor
            if qual > min_qual
                ϕᵢᴸ, ϕᵢᵁ = restrict(ϕ[i], U, k)
                sᵢᴸ, sᵢᵁ = sgn(sign(lb), s[i], false, -1), sgn(sign(lb), s[i], false, 1)
                push!(ϕ̃, ϕᵢᴸ, ϕᵢᵁ)
                push!(s̃, sᵢᴸ, sᵢᵁ)
            else
                # split along largest direction
                @debug "Splitting $U along $k"
                h, k = findmax(xu - xl)
                hmin = 1e-4 # FIXME: make a parameter?
                if h < hmin # stop splitting if the box is too small
                    @warn "Terminal case of recursion reached on $U, resorting to low-order method"
                    return f(xc) * prod(xu - xl)
                else
                    Uₗ, Uᵣ = split(U, k)
                    Ωₗ, Ωᵣ = ImplicitDomain(Uₗ, ϕ, s), ImplicitDomain(Uᵣ, ϕ, s)
                end
                return integrate(f, Ωₗ, quad) + integrate(f, Ωᵣ, quad)
            end
        end
        # If we got here then we have a good direction k do recurse down on the
        # dimension
        @debug "Recursing down on $k for $U"
        Ũ = remove_dimension(U, k)
        f̃ = integrand_eval(f, Ω, k, quad)
        Ω̃ = ImplicitDomain(Ũ, ϕ̃, s̃)
        return integrate(f̃, Ω̃, quad)
    end
end

"""
    sgn(m, s, S::Bool, σ)

Helper function to compute the sign of lower and upper restrictions of a
level-set function `ϕᵢ` in a box along a given height direction `k`. Here `m` is
sign of `∂ₖϕᵢ`, which is assume not to change throughout the box since `k` is
assumed to be a height direction, s is the sign of `ϕᵢ` in the box, S is a flag
to indicate whether we are in the special case of a surface integral.
"""
function sgn(m, s, S::Bool, σ)
    if (m == σ * s) || S
        return σ * m
    else
        return 0
    end
end

# export bound, restrict

end # module

"""
    struct Config

The `Config` struct represents the configuration for implicit integration,
passed to the [`integrate`](@ref) function to customize the behavior of the
algorithm.

The `Config` struct has the following fields:

- `quad`: called as `quad(f, a::SVector, b::SVector)`, returns an approximation
  to the integral of `f` over [`HyperRectangle(a,b)`(@ref)].
- `find_zero`: called as `find_zero(f, a, b)`, returns a zero of `f` in the
  interval `[a,b]` (if it exists).
- `find_zeros`: called as `find_zeros(f, a, b)`, returns all zeros of `f` in the
  interval `[a,b]` (if they exist).
- `min_qual`: a number between `0` and `1` used to specify the minimum quality
  factor for a height direction to be considered good for recursion. The quality
  factor for a direction `k` is given by `|∂ₖϕ| / |∇ϕ|`.
- `min_size`: a number used to specify the minimum size of a box to be split
  during the recursion. Recursion stops when the box is smaller than this size.
"""
@kwdef struct Config{T1,T2,T3}
    quad::T1          = (f, a, b) -> HCubature.hcubature(f, a, b)[1]
    find_zero::T2     = Roots.find_zero
    find_zeros::T3    = Roots.find_zeros
    min_qual::Float64 = 0.0
    min_size::Float64 = 1e-4
end

function integrate(
    f,
    ϕ,
    lc::SVector{N,T},
    hc::SVector{N,T};
    surface = false,
    config = Config(),
) where {N,T}
    U = HyperRectangle(lc, hc)
    RET_TYPE = typeof(f(lc)) # a guess for the return type...
    ϕ_ = SubFunction{N}(ϕ, SVector{0,Int}(), SVector{0,T}())
    s = surface ? 0 : -1
    return _integrate(f, [ϕ_], [s], U, config, RET_TYPE, surface)
end

function _integrate(
    f,
    phi_vec,
    s_vec,
    U::HyperRectangle{DIM,T},
    config,
    ::Type{RET_TYPE},
    surface::Bool,
) where {DIM,T,RET_TYPE}
    xl, xu = bounds(U)
    # Start by prunning phi_vec...
    partial_cell_idxs = Int[]
    for i in eachindex(phi_vec, s_vec)
        c = cell_type(phi_vec[i], s_vec[i], U, surface)
        c == empty_cell && return zero(RET_TYPE)
        c == partial_cell && push!(partial_cell_idxs, i)
    end
    if length(partial_cell_idxs) == 0
        return config.quad(f, xl, xu) # full cell
    end
    phi_vec = phi_vec[partial_cell_idxs]
    s_vec   = s_vec[partial_cell_idxs]
    # Finished prunning. If we did not return before this point, then the domain
    # is neither empty nor full. Next try to find a good direction to recurse
    # on. We will choose the direction with the largest gradient.
    if DIM == 1 # base case
        f̃ = _integrand_eval(f, phi_vec, s_vec, U, 1, config, RET_TYPE)
        x̃ = SVector{0,T}() # zero-argument vector to evaluate `f̃` (a const.)
        return f̃(x̃)
    end
    xc = (xl + xu) / 2
    ∇ϕ₁ = (x) -> gradient(phi_vec[1], x)
    k = argmax(abs.(∇ϕ₁(xc)))
    # Now check if k is a "good" height direction for all the level-set functions
    s_vec_new = Int[]
    R = restriction_type(eltype(phi_vec))
    phi_vec_new = R[]
    for i in eachindex(phi_vec, s_vec)
        ∇ϕᵢ_bnds = bound_gradient(phi_vec[i], U)
        den = sum(∇ϕᵢ_bnds) do (lb, ub)
            return max(abs(lb), abs(ub))^2
        end |> sqrt # max over U of |∇ϕᵢ|
        lb, ub = ∇ϕᵢ_bnds[k]
        qual = den == 0 ? 1.0 : lb * ub > 0 ? min(abs(lb), abs(ub)) / den : 0.0 # |∂ₖϕᵢ| / |∇ϕᵢ|
        if qual > config.min_qual
            # Restrict the level-set function to the box and push it to new list
            ϕᵢᴸ, ϕᵢᵁ = restrict(phi_vec[i], U, k)
            sign_∂ₖ = lb < 0 ? -1 : 1 # sign of ∂ₖϕᵢ
            sᵢᴸ, sᵢᵁ = sgn(sign_∂ₖ, s_vec[i], false, -1), sgn(sign_∂ₖ, s_vec[i], false, 1)
            push!(phi_vec_new, ϕᵢᴸ, ϕᵢᵁ)
            push!(s_vec_new, sᵢᴸ, sᵢᵁ)
        else
            # Direction k now good for recursion on dimension, so immediately
            # recurse on the box size
            h, dir = findmax(xu - xl)
            # split along largest direction
            @debug "Splitting $U along $dir"
            if h < config.min_size # stop splitting if the box is too small
                @warn "Terminal case of recursion reached on $U, resorting to low-order method."
                @debug "Tried to recurse along direction $k, but got a quality factor of $qual."
                return f(xc) * prod(xu - xl)
            else
                Uₗ, Uᵣ = split(U, dir)
            end
            return _integrate(f, phi_vec, s_vec, Uₗ, config, RET_TYPE, surface) +
                   _integrate(f, phi_vec, s_vec, Uᵣ, config, RET_TYPE, surface)
        end
    end
    # k is a good height direction for all the level-set functions, so recurse
    # on dimension until 1D integrals are reached
    @debug "Recursing down on $k for $U"
    if surface
        @assert length(phi_vec) == 1
        f̃ = _surface_integrand_eval(f, phi_vec[1], U, k, config, RET_TYPE)
    else
        f̃ = _integrand_eval(f, phi_vec, s_vec, U, k, config, RET_TYPE)
    end
    Ũ = remove_dimension(U, k)
    return _integrate(f̃, phi_vec_new, s_vec_new, Ũ, config, RET_TYPE, false)
end

## Algorithm 1 of Saye 2015
"""
    _integrand_eval(f, V::ImplicitDomain{N}, dir, quad1d)

Return a function `f̃ : ℝᴺ⁻¹ -> ℝ` that approximates the one-dimensional
integral over `I(x̃) ⊂ ℝ` of the function `t -> f(insert(x̃, k, t))`, where the
integration domain `I` is defined as `I(x̃) = {t ∈ [a,b] :
sᵢ*ϕᵢ(insert(̃x,k,t) ≥ 0 ∀ (ϕᵢ,sᵢ) ∈ V)}`.
"""
function _integrand_eval(
    f,
    phi_vec,
    s_vec,
    U::HyperRectangle{N},
    k::Int,
    config,
    ::Type{RET_TYPE},
) where {N,RET_TYPE}
    xl, xu = bounds(U)
    a, b = xl[k], xu[k]
    f̃ = (x̃) -> begin
        # compute the connected components
        bnds = [a, b]
        for ϕᵢ in phi_vec
            g = (t) -> ϕᵢ(insert(x̃, k, t))
            if N == 1
                # possible several zeros
                append!(bnds, config.find_zeros(g, (a, b)))
            else
                # we know that g is monotonic since it corresponds to a
                # height-direction, so at most a single root exists.
                # TODO: look up and test Brent's method for this case
                g(a) * g(b) > 0 && continue
                push!(bnds, config.find_zero(g, (a, b)))
            end
        end
        sort!(bnds)
        # compute the integral
        acc = zero(RET_TYPE)
        for i in 1:(length(bnds)-1) # loop over each segment
            rᵢ, rᵢ₊₁ = bnds[i], bnds[i+1]
            L = rᵢ₊₁ - rᵢ
            # decide if the segment is inside the domain
            xc = insert(x̃, k, (rᵢ + rᵢ₊₁) / 2)
            skip = false
            for (ϕᵢ, sᵢ) in zip(phi_vec, s_vec)
                sᵢ == 0 && continue # avoid evaluation of ϕᵢ(xc) when possible
                if sᵢ * ϕᵢ(xc) < 0
                    skip = true
                    break
                end
            end
            skip && continue
            # add the contribution of the segment by performing a 1D quadrature.
            acc += config.quad(SVector(rᵢ), SVector(rᵢ₊₁)) do (t,)
                x = insert(x̃, k, t)
                return f(x)
            end
        end
        return acc
    end
    return f̃
end

function _surface_integrand_eval(
    f,
    phi,
    U::HyperRectangle{N},
    k::Int,
    config,
    ::Type{RET_TYPE},
) where {N,RET_TYPE}
    xl, xu = bounds(U)
    a, b = xl[k], xu[k]
    f̃ = (x̃) -> begin
        g = (t) -> phi(insert(x̃, k, t))
        if g(a) * g(b) > 0
            return zero(RET_TYPE)
        else
            root = config.find_zero(g, (a, b))
            x = insert(x̃, k, root)
            ∇ϕ = gradient(phi, x)
            return f(x) * norm(∇ϕ) / abs(∇ϕ[k])
        end
    end
    return f̃
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
    v = σ * m
    if (m == σ * s) || S
        return v
    else
        return zero(v)
    end
end

"""
    @enum CellType

Enumeration for different types of cells. Options are `full_cell`, `empty_cell`,
and `partial_cell`.
"""
@enum CellType full_cell empty_cell partial_cell

"""
    cell_type(ϕ, s, U, surface)

Compute the [`CellType`](@ref) of a cell defined by the level-set function `ϕ`,
a sign `s`, and the box `U`. If `surface` is `true`, then the cell is classified
as per a surface integral.
"""
function cell_type(ϕ, s, U, surface)
    lb, ub = bound(ϕ, U)
    if lb * ub ≤ 0 # sign change
        return partial_cell
    elseif surface # no sign change, so no zero for surface integral
        return empty_cell
    else # no sign change, but since this is a volume integral cell can be full or empty
        if (lb > 0 && s < 0) || (ub < 0 && s > 0)
            return empty_cell
        else
            return full_cell
        end
    end
end

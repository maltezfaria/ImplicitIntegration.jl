"""
    struct Config

The `Config` struct represents the configuration for implicit integration, passed to the
[`integrate`](@ref) function to customize the behavior of the algorithm. It contains the
following fields:

  - `find_zero` a function with signature `(f, a, b, tol) --> x` such that `f(x) ≈ 0`, `a ≤ x ≤ b`. The tolerance `tol` is used to specify the absolute tolerance of the zero
    approximation (e.g. `xatol` in `Roots`).
  - `quad`: a function with signature `quad(f, a, b, tol) --> (I,E)` such that `I`
    approximates the integral of `f` over `[a,b]` and `E` is the estimated error. `a` and `b`
    `Tuple`(s)/`SVector`(s) specifying the lower and upper bounds of the integration domain,
    and `tol` is the desired absolute tolerance.
  - `min_vol`: a function with signature `(tol) --> Float64` used to specify the volume of
    boxes below which the spatial subdivision stops.
    low-order method is used to approximate the integral.
  - `min_qual`: a number between `0` and `1` used to specify the minimum quality factor for a
    height direction to be considered valid for recursion. The quality factor for a direction
    `k` is given by `|∂ₖϕ| / |∇ϕ|`.
"""
@kwdef struct Config{T1,T2,T3}
    find_zero::T1 = (f, a, b, tol) -> Roots.find_zero(f, (a, b), Roots.Brent(); xatol = tol)
    quad::T2 = (f, a, b, tol) -> HCubature.hcubature(f, a, b; atol = tol)
    min_vol::T3 = (tol) -> sqrt(eps(Float64))
    min_qual::Float64 = 0.1
end

"""
    struct LogInfo

A structure to store logging information for integration processes.

# Fields

  - `subdivisions::Vector{Int}`: A vector containing the subdivisions per
    dimension used during the integration process.
  - `loworder::Int`: The number of times the low-order method was used.
"""
mutable struct LogInfo
    subdivisions::Vector{Int}
    loworder::Int
    fullcells::Int
    tree::TreeNode
end

"""
    LogInfo(U::HyperRectangle)

Initialize a `LogInfo` object for integrating over `U`.
"""
function LogInfo(U::HyperRectangle{N}) where {N}
    subdivs = zeros(Int, N)
    loworder = 0
    fullcells = 0
    tree = TreeNode(U)
    return LogInfo(subdivs, loworder, fullcells, tree)
end

## overload the display method to better visualize the info
function Base.show(io::IO, info::LogInfo)
    println(io, "LogInfo:")
    println(io, "|- Subdivisions:")
    for (i, s) in enumerate(info.subdivisions)
        println(io, "|-- Dimension $i: $s subdivisions")
    end
    info.loworder > 0 && println(io, "|- Low-order used $(info.loworder) times")
    return io
end

"""
    integrate(f, ϕ, lc, hc; tol=1e-8, surface=false, config = Config(), loginfo = false) -->
    (val, logger)

Integrate the function `f` over an implicit domain defined by:

  - `Ω = {lc ≤ 𝐱 ≤ hc: ϕ(𝐱) < 0}` if `surface = false`
  - `Γ = {lc ≤ 𝐱 ≤ hc: ϕ(𝐱) = 0}` if `surface = true`

where `lc` and `hc` denote the lower and upper corners of the bounding box.

`tol` specifies the desired (absolute) tolerance for the approximation.

The function returns a named tuple `(val, logger)` where `val` is the approximated value,
and `logger` is a [`LogInfo`](@ref) object containing information about the integration
process if `loginfo = true`; otherwise, `logger` is `nothing`.

For a finer control over the integration process, pass a `config` object (see
[`Config`](@ref)).

Note that both `f` and `ϕ` must be callable with a single argument `𝐱` of type `SVector`.
Furthemore, `ϕ` is expected to return a real value.

See also [`quadgen`](@ref) if you want to generate a quadrature instead of directly computing
the value of the integral.

By default, `ImplicitIntegration` uses `ForwardDiff` to compute gradients and
`IntervalArithmetic` to compute bounds, both of which are needed for the algorithm to work.
While these work reasonably well in most cases, you may want to overload the following
methods for the type of your input function `ϕ`:

  - `ϕ(x::SVector{N,<:Real}) -> Real` to evaluate the level-set function at `x`.
  - `ϕ(xI::SVector{N,<:Interval{<:Real}}) -> Interval{<:Real}` to evaluate a bound on `ϕ` on
    the interval `xI`.
  - `ϕ(xD::SVector{N,Dual{N,<:Real}}) -> Dual{N,<:Real}` to evaluate the level-set function
    and its gradient at `x`.
  - `ϕ(xDI::SVector{N,Dual{N,<:Interval{<:Real}}}) -> Dual{N,<:Interval{<:Real}}` to
    evaluate a bound on `ϕ` and its gradient on the interval `xDI`.

You may need to overload the methods above if `typeof(ϕ)` is not supported by `ForwardDiff`
and/or `IntervalArithmetic`, or if you have a better/faster implementation; see the main docs.

# Examples

To compute the area of a quarter of a disk of radius 1.0:

```jldoctest; output = false
a, b = (0.0, 0.0), (1.5, 1.5)
ϕ = (x) -> x[1]^2 + x[2]^2 - 1
f = (x) -> 1.0
res = integrate(f, ϕ, a, b) # area of quarter of a disk
res.val ≈ π / 4

# output

true

```

To compute the perimeter of a quarter of a circle of radius 1.0:

```jldoctest; output = false
a, b = (0.0, 0.0), (1.5, 1.5)
ϕ = (x) -> x[1]^2 + x[2]^2 - 1
f = (x) -> 1.0
res = integrate(x -> 1.0, ϕ, a, b; surface = true) # perimeter of quarter of a circle
res.val ≈ 2π / 4

# output

true

```
"""
function integrate(
    f,
    ϕ,
    lc::SVector{N,T},
    hc::SVector{N,T};
    surface = false,
    tol = 1e-8,
    config = Config(),
    loginfo = false,
) where {N,T}
    U = HyperRectangle(lc, hc)
    logger = loginfo ? LogInfo(U) : nothing
    tree = loginfo ? logger.tree : nothing
    RET_TYPE = typeof(f(lc) * one(T) + f(hc) * one(T)) # a guess for the return type...
    s = surface ? 0 : -1
    val = _integrate(f, [ϕ], [s], U, config, RET_TYPE, Val(surface), tol, logger, tree)
    return (; val, logger)
end

function integrate(f, ϕ, lc, hc; kwargs...)
    @assert length(lc) == length(hc) "Lower and upper corners must have the same length."
    N = length(lc)
    T1, T2 = eltype(lc), eltype(hc)
    T = promote_type(float(T1), float(T2))
    return integrate(f, ϕ, SVector{N,T}(lc), SVector{N,T}(hc); kwargs...)
end

@noinline function _integrate(
    f,
    phi_vec,
    s_vec,
    U::HyperRectangle{DIM,T},
    config,
    ::Type{RTYPE},
    ::Val{S},
    tol,
    logger,
    tree,
) where {DIM,T,RTYPE,S}
    xl, xu = bounds(U)
    # Start by pruning phi_vec...
    partial_cell_idxs = Int[]
    for i in eachindex(phi_vec, s_vec)
        c = cell_type(phi_vec[i], s_vec[i], U, S)
        c == empty_cell && return zero(RTYPE)
        c == partial_cell && push!(partial_cell_idxs, i)
    end
    if length(partial_cell_idxs) == 0 # full cell
        isnothing(logger) || (logger.fullcells += 1)
        val, _ = config.quad(f, xl, xu, tol)
        return val
    end
    phi_vec = phi_vec[partial_cell_idxs]
    s_vec = s_vec[partial_cell_idxs]
    grad_phi_vec = map(gradient, phi_vec)
    # Finished pruning. If we did not return before this point, then the domain is neither
    # empty nor full. Next try to find a good direction to recurse on. We will choose the
    # direction with the largest gradient.
    if DIM == 1 # base case
        f̃ = if S
            @assert length(phi_vec) == 1
            _surface_integrand_eval(
                f,
                phi_vec[1],
                grad_phi_vec[1],
                U,
                1,
                config,
                RTYPE,
                tol,
                logger,
                tree,
            )
        else
            _integrand_eval(
                f,
                phi_vec,
                grad_phi_vec,
                s_vec,
                U,
                1,
                config,
                RTYPE,
                tol,
                logger,
                tree,
            )
        end
        x̃ = SVector{0,T}() # zero-argument vector to evaluate `f̃` (a const.)
        @debug "Reached 1D base case, evaluating integrand at $x̃"
        return f̃(x̃)
    end
    xc = (xl + xu) / 2
    ∇ϕ₁ = grad_phi_vec[1]
    k = argmax(abs.(∇ϕ₁(xc)))
    # Now check if k is a "good" height direction for all the level-set functions
    s_vec_new = Int[]
    R = Any # type of restriction. TODO: infer the type?
    phi_vec_new = R[]
    for i in eachindex(phi_vec, s_vec)
        ∇ϕᵢ_bnds = bound(grad_phi_vec[i], U)
        den = sum(∇ϕᵢ_bnds) do (lb, ub)
            return max(abs(lb), abs(ub))^2
        end |> sqrt # max over U of |∇ϕᵢ|
        lb, ub = ∇ϕᵢ_bnds[k]
        qual = den == 0 ? 1.0 : lb * ub > 0 ? min(abs(lb), abs(ub)) / den : 0.0 # |∂ₖϕᵢ| / |∇ϕᵢ|
        if qual > config.min_qual
            # Restrict the level-set function to the box and push it to new list
            ϕᵢᴸ, ϕᵢᵁ = project(phi_vec[i], k, xl[k]), project(phi_vec[i], k, xu[k])
            sign_∂ₖ = lb < 0 ? -1 : 1 # sign of ∂ₖϕᵢ
            sᵢᴸ, sᵢᵁ = sgn(sign_∂ₖ, s_vec[i], false, -1), sgn(sign_∂ₖ, s_vec[i], false, 1)
            push!(phi_vec_new, ϕᵢᴸ, ϕᵢᵁ)
            push!(s_vec_new, sᵢᴸ, sᵢᵁ)
        else
            # Direction k not good for recursion on dimension, so immediately
            # recurse on the box size
            _, dir = findmax(xu - xl)
            vol = S ? sum(xu - xl) : prod(xu - xl)
            # split along largest direction
            if vol < config.min_vol(tol) # stop splitting if the box is too small
                isnothing(logger) || (logger.loworder += 1)
                @warn "Terminal case of recursion reached on $U, resorting to low-order method."
                if !S && all(i -> phi_vec[i](xc) * s_vec[i] > 0, 1:length(phi_vec))
                    return f(xc) * prod(xu - xl)
                else
                    return zero(RTYPE)
                end
            else # split the box
                @debug "Splitting $U along $dir"
                Uₗ, Uᵣ = split(U, dir)
                # compute the restriction of the level-sets on the left and right
                phi_vec_left = empty(phi_vec)
                phi_vec_right = empty(phi_vec)
                for ϕ in phi_vec
                    ϕₗ, ϕᵣ = split(ϕ, U, dir)
                    push!(phi_vec_left, ϕₗ)
                    push!(phi_vec_right, ϕᵣ)
                end
                tree_left = isnothing(tree) ? nothing : TreeNode(Uₗ)
                tree_right = isnothing(tree) ? nothing : TreeNode(Uᵣ)
                isnothing(tree) || (push!(tree.children, (tree_left, 0), (tree_right, 0)))
                isnothing(logger) || (logger.subdivisions[DIM] += 1)
                tol /= 2 # FIXME: halving the tolerance is way too much in practice...
                Iₗ = _integrate(
                    f,
                    phi_vec_left,
                    s_vec,
                    Uₗ,
                    config,
                    RTYPE,
                    Val(S),
                    tol,
                    logger,
                    tree_left,
                )
                Iᵣ = _integrate(
                    f,
                    phi_vec_right,
                    s_vec,
                    Uᵣ,
                    config,
                    RTYPE,
                    Val(S),
                    tol,
                    logger,
                    tree_right,
                )
                return Iₗ + Iᵣ
            end
        end
    end
    # k is a good height direction for all the level-set functions, so recurse
    # on dimension until 1D integrals are reached
    @debug "Recursing down on $k for $U"
    Ũ = remove_dimension(U, k)
    subtree = isnothing(tree) ? nothing : TreeNode(Ũ)
    isnothing(tree) || (push!(tree.children, (subtree, k)))
    if S
        @assert length(phi_vec) == 1
        f̃ = _surface_integrand_eval(
            f,
            phi_vec[1],
            grad_phi_vec[1],
            U,
            k,
            config,
            RTYPE,
            tol,
            logger,
            tree,
        )
        return _integrate(
            f̃,
            phi_vec_new,
            s_vec_new,
            Ũ,
            config,
            RTYPE,
            Val(false),
            tol,
            logger,
            subtree,
        )
    else
        f̃ = _integrand_eval(
            f,
            phi_vec,
            grad_phi_vec,
            s_vec,
            U,
            k,
            config,
            RTYPE,
            tol,
            logger,
            tree,
        )
        return _integrate(
            f̃,
            phi_vec_new,
            s_vec_new,
            Ũ,
            config,
            RTYPE,
            Val(false),
            tol,
            logger,
            subtree,
        )
    end
end

## Algorithm 1 of Saye 2015
"""
    _integrand_eval(f, phi_vec, s_vec, U, k, config, ::Type{RET_TYPE}, tol, logger)

Return a function `f̃ : ℝᴺ⁻¹ -> ℝ` that approximates the one-dimensional integral over
`I(x̃) ⊂ ℝ` of the function `t -> f(insert(x̃, k, t))`, where the integration domain `I` is
defined as `I(x̃) = {t ∈ [a,b] : sᵢ*ϕᵢ(insert(̃x,k,t) ≥ 0 ∀ (ϕᵢ,sᵢ) ∈ V)}`.
"""
function _integrand_eval(
    f,
    phi_vec,
    grad_phi_vec,
    s_vec,
    U::HyperRectangle{N},
    k::Int,
    config,
    ::Type{RET_TYPE},
    tol,
    logger,
    tree,
) where {N,RET_TYPE}
    xl, xu = bounds(U)
    a, b = xl[k], xu[k]
    f̃ = (x̃) -> begin
        # compute the connected components
        bnds = [a, b]
        for (ϕᵢ, ∇ϕᵢ) in zip(phi_vec, grad_phi_vec)
            if N == 1
                # possible several zeros. Use internal `find_zeros` method which
                # works on the function `ϕᵢ` directly so that it can tap into
                # the `bound(ϕᵢ)` and `bound(∇ϕᵢ)` methods.
                _find_zeros!(bnds, ϕᵢ, ∇ϕᵢ, U, config, tol, logger, tree)
            else
                # we know that g is monotonic since it corresponds to a
                # height-direction, so at most a single root exists.
                g = (t) -> ϕᵢ(insert(x̃, k, t))
                g(a) * g(b) > 0 && continue
                push!(bnds, config.find_zero(g, a, b, tol))
            end
        end
        sort!(bnds)
        # HACK: keep only the unique elements in bnds up to ≈ 1e-8 relative tolerance.
        # Avoids cases where we have two zeros that are very close to each other, usually
        # coming from a degenerate root (e.g. x^2 at x = 0 with numerical noise).
        unique!(x -> round(x; sigdigits = 8), bnds)
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
            val, _ = config.quad(SVector(rᵢ), SVector(rᵢ₊₁), tol) do (t,)
                x = insert(x̃, k, t)
                return f(x)
            end
            acc += val
        end
        return acc
    end
    return f̃
end

function _surface_integrand_eval(
    f,
    phi,
    phi_grad,
    U::HyperRectangle{N,T},
    k::Int,
    config,
    ::Type{RET_TYPE},
    tol,
    logger,
    tree,
) where {N,T,RET_TYPE}
    xl, xu = bounds(U)
    a, b = xl[k], xu[k]
    f̃ = (x̃) -> begin
        g = (t) -> phi(insert(x̃, k, t))
        if N == 1
            # corner case where we have a "surface" integral in 1D. Arises only when calling
            # `integrate` with `surface=true` on one-dimensional level-set functions.
            roots = T[]
            _find_zeros!(roots, phi, phi_grad, U, config, tol, logger, tree)
            sum(roots) do root
                x = insert(x̃, k, root)
                ∇ϕ = phi_grad(x)
                return f(x) * norm(∇ϕ) * inv(abs(∇ϕ[k]))
            end
        else # guaranteed to have at most one zero
            if g(a) * g(b) > 0
                return zero(RET_TYPE)
            else
                root = config.find_zero(g, a, b, tol)
                x = insert(x̃, k, root)
                ∇ϕ = phi_grad(x)
                return f(x) * norm(∇ϕ) * inv(abs(∇ϕ[k]))
            end
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

"""
    _find_zeros!(roots, f, U::Segment, config, tol)

Return all zeros of the function `f` in the `Segment` `U`. `f` should be callable as
`f(x::SVector{1})`.
"""
function _find_zeros!(roots, ϕ, ∇ϕ, U::Segment, config, tol, logger, tree)
    xl, xu = bounds(U)
    if norm(xu - xl) < config.min_vol(tol)
        # splitting has led to very small boxes, likely due to e.g. degenerate
        # roots (e.g. x^2 on [-1,1]). Give up on trying to find all zeros and
        # simply try a bisection.
        if ϕ(xl) * ϕ(xu) ≤ 0
            g = (t) -> ϕ(SVector(t))
            r = config.find_zero(g, xl[1], xu[1], tol)
            push!(roots, r)
            return roots
        else # no zeros
            return roots
        end
    end

    ϕl, ϕu = bound(ϕ, U)
    if ϕl * ϕu > 0 # no zeros in the interval
        return roots
    else # maybe there are zeros
        ∇ϕl, ∇ϕu = bound(∇ϕ, U) |> first
        if ∇ϕl * ∇ϕu > 0 # monotonic, so at most one zero
            if ϕ(xl) * ϕ(xu) ≤ 0
                g = (t) -> ϕ(SVector(t))
                r = config.find_zero(g, xl[1], xu[1], tol)
                push!(roots, r)
                return roots
            else
                return roots
            end
        else # can't prove monotonicity nor lack of zeros, so split
            Uₗ, Uᵣ = split(U, 1)
            ϕₗ, ϕᵣ = split(ϕ, U, 1)
            tree_left = isnothing(tree) ? nothing : TreeNode(Uₗ)
            tree_right = isnothing(tree) ? nothing : TreeNode(Uᵣ)
            isnothing(tree) || (push!(tree.children, (tree_left, 0), (tree_right, 0)))
            isnothing(logger) || (logger.subdivisions[1] += 1) # one-dimensional subdivision
            _find_zeros!(roots, ϕₗ, gradient(ϕₗ), Uₗ, config, tol, logger, tree)
            _find_zeros!(roots, ϕᵣ, gradient(ϕᵣ), Uᵣ, config, tol, logger, tree)
            return roots
        end
    end
end

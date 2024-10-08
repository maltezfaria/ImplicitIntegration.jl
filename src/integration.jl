"""
    struct Config

The `Config` struct represents the configuration for implicit integration,
passed to the [`integrate`](@ref) function to customize the behavior of the
algorithm.

The `Config` struct has the following fields:

- `find_zero`: called as `find_zero(f, a, b)`, returns a zero of `f` in the
  interval `[a,b]` (if it exists).
- `quad`: called as `quad(f, a::SVector, b::SVector)`, returns an approximation
  to the integral of `f` over [`HyperRectangle(a,b)`(@ref)].
- `min_qual`: a number between `0` and `1` used to specify the minimum quality
  factor for a height direction to be considered good for recursion. The quality
  factor for a direction `k` is given by `|∂ₖϕ| / |∇ϕ|`.
- `min_vol`: a number used to specify the minimum size of a box to be split
  during the recursion. Recursion stops when the box's volume is smaller than
  this size. For surface integrals, the size is taken to be the area of the box.
"""
@kwdef struct Config{T1,T2}
    find_zero::T1     = (f, interval) -> Roots.find_zero(f, interval, Roots.Brent())
    quad::T2          = nothing
    min_qual::Float64 = 0.0
    min_vol::Float64  = 1e-8
end

"""
    Config(tol)

Create a `[Config](@ref)` object where all fields try to respect the given
tolerance `tol`. In practice, the following heuristics are used:

- `find_zero` uses `xatol=tol` so the implicit surface `ϕ = 0` is only sought to
  within `tol` of the true surface
- `min_vol` is set to `tol` so that the recursion stops when the integration
  measure of the box is smaller than `tol`
"""
function Config(tol::Real)
    return Config(;
        find_zero = (f, interval) ->
            Roots.find_zero(f, interval, Roots.Brent(); xatol = tol),
        min_vol = tol,
    )
end

"""
    struct LogInfo

A structure to store logging information for integration processes.

# Fields
- `subdivisions::Vector{Int}`: A vector containing the subdivisions per
  dimension used during the integration process.
"""
struct LogInfo
    subdivisions::Vector{Int}
end

"""
    LogInfo(d::Integer)

Initialize a `LogInfo` object for integrating over `ℝᵈ`.
"""
function LogInfo(dim::Integer)
    subdivs = zeros(Int, dim)
    return LogInfo(subdivs)
end

## overload the display  method to better visualize the info
function Base.show(io::IO, ::MIME"text/plain", info::LogInfo)
    println(io, "LogInfo:")
    println(io, "|- Subdivisions:")
    for (i, s) in enumerate(info.subdivisions)
        println(io, "|-- Dimension $i: $s subdivisions")
    end
end

"""
    integrate(f, ϕ, lc, hc; tol=1e-8, surface=false, config = Config())

Integrate the function `f` over an implict domain defined by:

- `Ω = {lc ≤ 𝐱 ≤ hc: ϕ(𝐱) < 0}` if `surface = false`
- `Γ = {lc ≤ 𝐱 ≤ hc: ϕ(𝐱) = 0}` if `surface = true`

where `lc::NTuple` and `hc::NTuple` denote the lower and upper corners of the bounding box.

`tol` specifies the desired (absolute) tolerance for the approximation.

For a finer control of the integration process, the user can pass a `config`
object to customize the behavior of various aspects of the algorithm (see
[`Config`](@ref) for more details). In such cases the `tol` parameter is ignored
and `config.quad` is used for the integration.

Note that both `f` and `ϕ` must be callable with a single argument `𝐱` of type
`SVector`. Furthemore, `ϕ` is expected to return a real value.

See also [`quadgen`](@ref) if you want to generate a quadrature instead of
direcly computing the value of the integral.

# Examples

To compute the area of a quarter of a disk of radius 1.0:

```jldoctest; output = false
a, b = (0.0, 0.0), (1.5, 1.5)
ϕ = (x) -> x[1]^2 + x[2]^2 - 1
f = (x) -> 1.0
integrate(f, ϕ, a, b) ≈ π / 4 # area of quarter of a disk

# output

true

```

To compute the perimeter of a quarter of a circle of radius 1.0:

```jldoctest; output = false
a, b = (0.0, 0.0), (1.5, 1.5)
ϕ = (x) -> x[1]^2 + x[2]^2 - 1
f = (x) -> 1.0
integrate(x -> 1.0, ϕ, a, b; surface = true) ≈ 2π / 4 # perimeter of quarter of a circle

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
    config = Config(tol),
    log = false,
) where {N,T}
    U = HyperRectangle(lc, hc)
    RET_TYPE = typeof(f(lc) * one(T) + f(hc) * one(T)) # a guess for the return type...
    ϕ_ = SubFunction{N}(ϕ, SVector{0,Int}(), SVector{0,T}())
    s = surface ? 0 : -1
    logger = log ? LogInfo(N) : nothing
    val = _integrate(f, [ϕ_], [s], U, config, RET_TYPE, Val(surface), tol, logger)
    return log ? (val, logger) : val
end

function integrate(f, ϕ, lc, hc; kwargs...)
    @assert length(lc) == length(hc) "Lower and upper corners must have the same length."
    N = length(lc)
    T1, T2 = eltype(lc), eltype(hc)
    T = promote_type(float(T1), float(T2))
    return integrate(f, ϕ, SVector{N,T}(lc), SVector{N,T}(hc); kwargs...)
end

function _integrate(
    f,
    phi_vec,
    s_vec,
    U::HyperRectangle{DIM,T},
    config,
    ::Type{RTYPE},
    ::Val{S},
    tol,
    logger,
) where {DIM,T,RTYPE,S}
    xl, xu = bounds(U)
    # Start by prunning phi_vec...
    partial_cell_idxs = Int[]
    for i in eachindex(phi_vec, s_vec)
        c = cell_type(phi_vec[i], s_vec[i], U, S)
        c == empty_cell && return zero(RTYPE)
        c == partial_cell && push!(partial_cell_idxs, i)
    end
    if length(partial_cell_idxs) == 0 # full cell
        if isnothing(config.quad)
            return HCubature.hcubature(f, xl, xu; atol = tol)[1]
        else
            return config.quad(f, xl, xu)
        end
    end
    phi_vec = phi_vec[partial_cell_idxs]
    s_vec   = s_vec[partial_cell_idxs]
    # Finished prunning. If we did not return before this point, then the domain
    # is neither empty nor full. Next try to find a good direction to recurse
    # on. We will choose the direction with the largest gradient.
    if DIM == 1 # base case
        f̃ = _integrand_eval(f, phi_vec, s_vec, U, 1, config, RTYPE, tol, logger)
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
            _, dir = findmax(xu - xl)
            vol    = S ? sum(xu - xl) : prod(xu - xl)
            # split along largest direction
            @debug "Splitting $U along $dir"
            if vol < config.min_vol # stop splitting if the box is too small
                @warn "Terminal case of recursion reached on $U, resorting to low-order method."
                @debug "Tried to recurse along direction $k, but got a quality factor of $qual."
                return f(xc) * prod(xu - xl)
            else
                Uₗ, Uᵣ = split(U, dir)
                isnothing(logger) || (logger.subdivisions[DIM] += 1)
            end
            tol′ = isnothing(tol) ? nothing : tol / 2
            return _integrate(f, phi_vec, s_vec, Uₗ, config, RTYPE, Val(S), tol′, logger) +
                   _integrate(f, phi_vec, s_vec, Uᵣ, config, RTYPE, Val(S), tol′, logger)
        end
    end
    # k is a good height direction for all the level-set functions, so recurse
    # on dimension until 1D integrals are reached
    @debug "Recursing down on $k for $U"
    Ũ = remove_dimension(U, k)
    if S
        @assert length(phi_vec) == 1
        f̃ = _surface_integrand_eval(f, phi_vec[1], U, k, config, RTYPE)
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
        )
    else
        f̃ = _integrand_eval(f, phi_vec, s_vec, U, k, config, RTYPE, tol, logger)
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
        )
    end
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
    tol,
    logger,
) where {N,RET_TYPE}
    xl, xu = bounds(U)
    a, b = xl[k], xu[k]
    quad = if isnothing(config.quad)
        (f, a, b) -> HCubature.hcubature(f, a, b; atol = tol)[1]
    else
        config.quad
    end
    f̃ = (x̃) -> begin
        # compute the connected components
        bnds = [a, b]
        for ϕᵢ in phi_vec
            if N == 1
                # possible several zeros. Use internal `find_zeros` method which
                # works on the function `ϕᵢ` directly so that it can tap into
                # the `bound` and `bound_gradient` methods.
                find_zeros!(bnds, ϕᵢ, U, config, logger)
            else
                # we know that g is monotonic since it corresponds to a
                # height-direction, so at most a single root exists.
                # TODO: look up and test Brent's method for this case
                g = (t) -> ϕᵢ(insert(x̃, k, t))
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
            acc += quad(SVector(rᵢ), SVector(rᵢ₊₁)) do (t,)
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
            return f(x) * norm(∇ϕ) * inv(abs(∇ϕ[k]))
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
    find_zeros(f, U::HyperRectangle{<:,T}, config = Config()) where {T}

Return all zeros of the function `f` in the box `U`. `f` should be callable as
`f(x::SVector)`.

The `config` argument can be used to customize the behavior of the algorithm by
passing e.g. a different `find_zero` method (e.g. `Roots.find_zero`) and a
different `min_vol` to control the box size at which recursion gives up trying
to perform a high-order integration and resorts to a midpoint rule.
"""
function find_zeros(
    f,
    U::HyperRectangle{<:,T},
    config = Config(),
    logger = nothing,
) where {T}
    return find_zeros!(T[], f, U, config, logger)
end

function find_zeros!(roots, ϕ, U::Segment, config, logger)
    xl, xu = bounds(U)
    if norm(xu - xl) < config.min_vol
        # splitting has led to very small boxes, likely due to e.g. degenerate
        # roots (e.g. x^2 on [-1,1]). Give up on trying to find all zeros and
        # simply try a bisection.
        if ϕ(xl) * ϕ(xu) ≤ 0
            g = (t) -> ϕ(SVector(t))
            r = config.find_zero(g, (xl[1], xu[1]))
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
        ∇ϕl, ∇ϕu = bound_gradient(ϕ, U) |> first
        if ∇ϕl * ∇ϕu > 0 # monotonic, so at most one zero
            if ϕ(xl) * ϕ(xu) ≤ 0
                g = (t) -> ϕ(SVector(t))
                r = config.find_zero(g, (xl[1], xu[1]))
                push!(roots, r)
                return roots
            else
                return roots
            end
        else # can't prove monotonicity nor lack of zeros, so split
            U1, U2 = split(U, 1)
            isnothing(logger) || (logger.subdivisions[1] += 1) # one-dimensional subdivision
            find_zeros!(roots, ϕ, U1, config, logger)
            find_zeros!(roots, ϕ, U2, config, logger)
            return roots
        end
    end
end

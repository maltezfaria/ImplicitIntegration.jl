struct GaussLegendre{N,T}
    nodes::SVector{N,T}
    weights::SVector{N,T}
end

"""
    GaussLegendre(; order, T = Float64)

Construct a Gauss-Legendre quadrature rule of given order with nodes and weights
of type `T`, callable as `quad1d(f, a, b)` where `f` is the function to
integrate and `a` and `b` are the bounds of the integration interval.

# Examples

```jldoctest; output = false
quad1d = ImplicitIntegration.GaussLegendre(;order = 20)
quad1d(cos,0,1) ≈ sin(1)

# output

true

```
"""
function GaussLegendre(; order, T::Type = Float64)
    N = ceil(Int, (order + 1) / 2)
    nodes_, weights_ = HCubature.QuadGK.gauss(T, N, 0, 1)
    return GaussLegendre(SVector{N,T}(nodes_), SVector{N,T}(weights_))
end

function (quad::GaussLegendre{N})(f, a::Number, b::Number) where {N}
    h = b - a
    acc = sum(zip(quad.nodes, quad.weights)) do (x̂, ŵ)
        x = a + h * x̂
        return f(x) * ŵ
    end
    return h * acc
end

"""
    TensorQuadrature(quad1d)

Given a 1D quadrature rule `quad1d` with signature `quad1d(f, a::Number,
b::Number)`, return a tensor product quadrature rule callable through
`quadnd(f,a::SVector{N}, b::SVector{N})`.
"""
struct TensorQuadrature{T}
    quad1d::T
end

function (Q::TensorQuadrature)(f, a::SVector{N}, b::SVector{N}) where {N}
    if N == 0 # base case
        return f(a)
    else
        Q(pop(a), pop(b)) do x̃
            Q.quad1d(a[N], b[N]) do xN
                x = push(x̃, xN)
                return f(x)
            end
        end
    end
end

"""
    struct Quadrature{N,T}

A collection of `coords` and `weights` for integration in `N` dimensions.

`Quadrature`s are the result of [`quadgen`](@ref) and can be used to integrate functions
through [`integrate(f,::Quadrature)`](@ref).
"""
struct Quadrature{N,T}
    coords::Vector{SVector{N,T}}
    weights::Vector{T}
end

"""
    integrate(f, Q::Quadrature)

Shorthand for `∑ᵢf(xᵢ)wᵢ`, where `xᵢ` and `wᵢ` are the nodes and weights of the
`Quadrature`.
"""
function integrate(f, Q::Quadrature{N,T}) where {N,T}
    if isempty(Q.coords)
        # return the zero of a plausible type when empty
        return f(zero(SVector{N,T})) * zero(T)
    else
        sum(zip(Q.coords, Q.weights)) do (x, w)
            return f(x) * w
        end
    end
end

# The code below is a hacky way to reuse the `integrate` to generate a
# quadrature instead of computing a sum. This is done by overloading the `*` and
# `+` operator to mean "multiply weight" and "concatenate quadratures"
# respectively. Not the most elegant solution, but it works. Note that the
# operations are in place to avoid unnecessary allocations, but this violates
# the semantics of e.g. `+` and `*`.

function Base.:(*)(w::Number, Q::Quadrature)
    rmul!(Q.weights, w)
    return Q
end
Base.:(*)(Q::Quadrature, w::Number) = w * Q

function Base.:(+)(Q::Quadrature, P::Quadrature)
    if length(P.coords) > length(Q.coords)
        Q, P = P, Q
    end
    append!(Q.coords, P.coords)
    append!(Q.weights, P.weights)
    return Q
end

Base.zero(::Type{Quadrature{N,T}}) where {N,T} = Quadrature(SVector{N,T}[], T[])
Base.zero(q::Quadrature) = zero(typeof(q))

"""
    quadgen(ϕ, lc, hc; order, surface=false, config = nothing)

Return a [`Quadrature`](@ref) to integrate a function over an implict domain
defined by:

- `Ω = {lc ≤ 𝐱 ≤ hc: ϕ(𝐱) < 0}` if `surface = false`
- `Γ = {lc ≤ 𝐱 ≤ hc: ϕ(𝐱) = 0}` if `surface = true`

where `lc::NTuple` and `hc::NTuple` denote the lower and upper corners of the bounding box.

The `order` parameter specifies the degree of exactness of the quadrature rule;
that is, the quadrature rule will integrate exactly polynomials of degree
up to `order`, but not `order+1`. A `GaussLegendre` quadrature rule is used.

For a finer control of the integration process, the user can pass a `config`
object to customize the behavior of various aspects of the algorithm (see
[`Config`](@ref) for more details). In such cases, the `order` parameter is
ignored and `config.quad` is used for the integration.

Note that `ϕ` must be callable with a single argument `𝐱` of type `SVector`.
Furthemore, `ϕ` is expected to return a real value.

# Examples

To compute the area of a quarter of a disk of radius 1.0:

```jldoctest; output = false
a, b = (0.0, 0.0), (1.5, 1.5)
ϕ = (x) -> x[1]^2 + x[2]^2 - 1
f = (x) -> 1.0
Q = quadgen(ϕ, a, b; order = 20)
integrate(f, Q) ≈ π / 4 # area of quarter of a disk

# output

true

```

To compute the perimeter of a quarter of a circle of radius 1.0:

```jldoctest; output = false
a, b = (0.0, 0.0), (1.5, 1.5)
ϕ = (x) -> x[1]^2 + x[2]^2 - 1
f = (x) -> 1.0
Q = quadgen(ϕ, a, b; order = 20, surface = true)
integrate(f, Q) ≈ 2π / 4 # perimeter of quarter of a circle

# output

true

```


"""
function quadgen(ϕ, lc::SVector{N,T}, hc::SVector{N,T}; order, kwargs...) where {N,T}
    if !haskey(kwargs, :config)
        quad1d = GaussLegendre(; order = order)
        quad = TensorQuadrature(quad1d)
        config = Config(; quad)
    else
        isnothing(order) || @warn "Ignoring `order` parameter since `config` is provided."
    end
    f = function (x::SVector{N,T}) where {N,T}
        # FIXME: this is possibly a peformance hog since it allocates new arrays
        # for each call. A more efficient way to handle could be to have instead
        # a `QNode` type that holds the current node and the current weight, and
        # overload operations on it.
        return Quadrature([x], [one(T)])
    end
    return integrate(f, ϕ, lc, hc; tol = nothing, config, kwargs...)
end

function quadgen(ϕ, lc, hc, args...; kwargs...)
    @assert length(lc) == length(hc) "Lower and upper corners must have the same length."
    N = length(lc)
    T1, T2 = eltype(lc), eltype(hc)
    T = promote_type(float(T1), float(T2))
    return quadgen(ϕ, SVector{N,T}(lc), SVector{N,T}(hc), args...; kwargs...)
end

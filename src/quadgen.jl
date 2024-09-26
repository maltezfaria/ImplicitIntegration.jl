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
qrule1d = GaussLegendre(;order = 5)
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

struct QNode{N,T}
    coord::SVector{N,T}
    weight::T
end

Base.:(*)(q::QNode, w::Number) = QNode(q.coord, q.weight * w)
Base.:(*)(w::Number, q::QNode) = q * w

"""
    struct Quadrature{N,T}

A collection of `coords` and `weights` for integration in `N` dimensions.

`Quadrature`s are the result of `quadgen` and can be used to integrate functions
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
function integrate(f, Q::Quadrature)
    @assert !isempty(Q.coords) "Quadrature must have at least one node."
    sum(zip(Q.coords, Q.weights)) do (x, w)
        return f(x) * w
    end
end

function Base.:(*)(w::Number, Q::Quadrature)
    rmul!(Q.weights, w)
    return Q
end
Base.:(*)(Q::Quadrature, w::Number) = w * Q

function Base.:(+)(Q::Quadrature, q::QNode)
    return (push!(Q.coords, q.coord); push!(Q.weights, q.weight); Q)
end
Base.:(+)(q::QNode, Q::Quadrature) = Q + q

function Base.:(+)(Q::Quadrature, P::Quadrature)
    if length(P.coords) > length(Q.coords)
        Q, P = P, Q
    end
    append!(Q.coords, P.coords)
    append!(Q.weights, P.weights)
    return Q
end

Base.:(+)(q::QNode, p::QNode) = Quadrature([q.coord, p.coord], [q.weight, p.weight])

Base.zero(::Type{Quadrature{N,T}}) where {N,T} = Quadrature(SVector{N,T}[], T[])
Base.zero(q::Quadrature) = zero(typeof(q))

function quadgen(ϕ, lc::SVector{N,T}, hc::SVector{N,T}; order, surface = false) where {N,T}
    quad1d = GaussLegendre(; order = order)
    quad   = TensorQuadrature(quad1d)
    config = Config(; quad)
    f      = function (x::SVector{N,T}) where {N,T}
        # return QNode(x, one(T))
        return Quadrature([x], [one(T)])
    end
    return integrate(f, ϕ, lc, hc; tol = nothing, surface, config)
end

function quadgen(ϕ, lc, hc; kwargs...)
    @assert length(lc) == length(hc) "Lower and upper corners must have the same length."
    N = length(lc)
    T1, T2 = eltype(lc), eltype(hc)
    T = promote_type(float(T1), float(T2))
    return quadgen(ϕ, SVector{N,T}(lc), SVector{N,T}(hc); kwargs...)
end


struct BernsteinPolynomial{N,T} <: Function
    coeffs::Array{T,N}
    low_corner::SVector{N,T}
    high_corner::SVector{N,T}
end

"""
    BernsteinPolynomial(c::AbstractArray, lc, hc)

Create a multidimensional [Bernstein polynomial](https://en.wikipedia.org/wiki/Bernstein_polynomial#Generalizations_to_higher_dimension) with coefficients `c` defined on the hyperrectangle
`[lc[1], hc[1]] × … × [lc[N], hc[N]]`.

Calling `p(x)` evaluates the polynomial at the point `x = (x[1], …, x[N])` by the formula

```math
p(x_1,\\dots,x_D)=\\sum_{i_j=0}^{d_j}c_{i_1\\dots i_D}\\prod_{j=1}^D\\binom{d_j}{i_j}(x_j-l_j)^{i_j}(r_j-x_j)^{d_j-i_j}
```

where ``l_j = lc[j]`` and ``r_j = hc[j]`` are the lower and upper bounds of the
hyperrectangle, respectively, and ``d_j = size(c)[j] - 1`` is the degree of the polynomial
in dimension `j`.

See also [`berninterp`](@ref).
"""
function BernsteinPolynomial(c::AbstractArray, lc, hc)
    N = ndims(c)
    return BernsteinPolynomial{N,eltype(c)}(c, SVector{N}(lc), SVector{N}(hc))
end

coefficients(p::BernsteinPolynomial) = p.coeffs
low_corner(p::BernsteinPolynomial)   = p.low_corner
high_corner(p::BernsteinPolynomial)  = p.high_corner
degree(p::BernsteinPolynomial)       = size(coefficients(p)) .- 1

# evaluation
function (p::BernsteinPolynomial{N})(x) where {N}
    x_ = SVector{N}(x) # try conversion to SVector
    l = low_corner(p)
    r = high_corner(p)
    x₀ = (x_ - l) ./ (r - l)
    c = coefficients(p)
    return _evaluate_bernstein(x₀, c, Val{N}(), 1, length(c))
end
@fastmath function _evaluate_bernstein(
    x::SVector{N},
    c::AbstractArray,
    ::Val{dim},
    i1,
    len,
) where {N,dim}
    n = size(c, dim)
    @inbounds xd = x[dim]
    # inspired by https://personal.math.ubc.ca/~cass/graphics/text/www/pdf/a6.pdf and the
    # FastChebInterp.jl package
    if dim == 1
        s = 1 - xd
        @inbounds P = c[i1]
        C = (n - 1) * xd
        for k in 1:(n-1)
            @inbounds P = P * s + C * c[i1+k]
            C = C * (n - k - 1) / (k + 1) * xd
        end
        return P
    else
        Δi = len ÷ n # column-major stride of current dimension

        # we recurse downward on dim for cache locality,
        # since earlier dimensions are contiguous
        dim′ = Val{dim - 1}()

        s = 1 - xd
        P = _evaluate_bernstein(x, c, dim′, i1, Δi)
        C = (n - 1) * xd
        for k in 1:(n-1)
            P = P * s + C * _evaluate_bernstein(x, c, dim′, i1 + k * Δi, Δi)
            C = C * (n - k - 1) / (k + 1) * xd
        end
        return P
    end
end

"""
    derivative(p::BernsteinPolynomial, d::Int)g

Compute the derivative along dimension `d` of the Bernstein polynomial `p`, returning a new
`BernsteinPolynomial` of the same dimension `N`.
"""
function derivative(p::BernsteinPolynomial{N}, d::Int) where {N}
    @assert 1 ≤ d ≤ N "Dimension $d out of bounds for polynomial of dimension $N"
    c = coefficients(p)
    n = size(c)[d]
    k = degree(p)[d]
    l = low_corner(p)
    u = high_corner(p)
    c′ = mapslices(c; dims = d) do b
        b′ = (b[2:n] .- b[1:(n-1)]) .* k ./ (u[d] - l[d])
        if k ≥ n
            push!(b′, -b[n] * k / (u[d] - l[d]))
        end
        return b′
    end
    return BernsteinPolynomial(c′, l, u)
end

"""
    split(p::BernsteinPolynomial, d::Integer, α = 0.5) --> pₗ, pᵣ

Split the Bernstein polynomial `p` along dimension `d` at `lc[d] + (hc[d] - lc[d]) * α`,
where `lc` and `hc` are the lower and upper corners of the hyperrectangle on which `p` is
defined. Returns two new Bernstein polynomials `pₗ` and `pᵣ` representing the left and right
polynomials.
"""
function split(p::BernsteinPolynomial{D,T}, d::Integer, α = 0.5) where {D,T}
    @assert 1 ≤ d ≤ D
    c = coefficients(p)
    k = degree(p)[d]
    k == 0 && return p, p
    n = size(c)[d]
    coeffs = mapslices(c; dims = d) do b
        c1 = Vector{T}()
        c2 = Vector{T}()
        for i in k:-1:1
            if i ≥ n
                push!(c1, b[1])
                @. b[1:(n-1)] = b[1:(n-1)] * (1 - α) + b[2:n] * α
                b[n] = b[n] * (1 - α)
            else
                push!(c1, b[1])
                pushfirst!(c2, b[i+1])
                @. b[1:i] = b[1:i] * (1 - α) + b[2:(i+1)] * α
            end
        end
        push!(c1, b[1])
        append!(c1, c2)
        return c1
    end
    lc, hc = low_corner(p), high_corner(p)
    split_pos = lc[d] + (hc[d] - lc[d]) * α
    p1 = BernsteinPolynomial(
        collect(selectdim(coeffs, d, 1:(k+1))),
        lc,
        setindex(hc, split_pos, d),
    )
    p2 = BernsteinPolynomial(
        collect(selectdim(coeffs, d, (k+1):(k+n))),
        setindex(lc, split_pos, d),
        hc,
    )
    return p1, p2
end

"""
    lower_restrict(p::BernsteinPolynomial{D}, d::Integer) where {D}

Restricts the given `BernsteinPolynomial` `p` to the lower face (i.e., where the `d`-th coordinate is at its lower bound)
along the specified dimension `d`. Returns a new `BernsteinPolynomial` of dimension `D-1` representing the restriction.
"""
function lower_restrict(p::BernsteinPolynomial{D,T}, d::Integer) where {D,T}
    @assert 1 ≤ d ≤ D
    c = coefficients(p)
    lc, hc = low_corner(p), high_corner(p)
    return BernsteinPolynomial(
        collect(selectdim(c, d, 1))::Array{T,D - 1},
        deleteat(lc, d),
        deleteat(hc, d),
    )
end

"""
    upper_restrict(p::BernsteinPolynomial{D}, d::Integer) where {D}

Same as [`lower_restrict`](@ref), but restricts to the upper face (i.e., where the `d`-th
coordinate is at its upper bound)
"""
function upper_restrict(p::BernsteinPolynomial{D,T}, d::Integer) where {D,T}
    @assert 1 ≤ d ≤ D
    c = coefficients(p)
    lc, hc = low_corner(p), high_corner(p)
    return BernsteinPolynomial(
        collect(selectdim(c, d, size(c)[d]))::Array{T,D - 1},
        deleteat(lc, d),
        deleteat(hc, d),
    )
end

"""
    bound(p::BernsteinPolynomial)

Return the bounds of the Bernstein polynomial `p` as a tuple `(m, M)`, where `m` is a lower
bound and `M` is an upper bound on the polynomial's values over the hyperrectangle with
corners `low_corner(p)` and `high_corner(p)`.
"""
function bound(p::BernsteinPolynomial)
    c = coefficients(p)
    m, M = extrema(c)
    c = coefficients(p)
    if prod(size(c)) < prod(degree(p) .+ 1)
        # some coefficients are implicitly zero, so take them into account in the bounds
        m = min(m, 0)
        M = max(M, 0)
    end
    return m, M
end

# interpolation utilities
"""
    uniform_points(n[, lc, hc])

Generate `n[1] × … × n[D]` uniformly spaced points in the
hyperrectangle defined by `(lc[1], hc[1]) × … × (lc[D], hc[D])`, where `D = length(n) = length(lc) = length(hc)`. By default `lc = (0, … , 0)` and `hc = (1, … , 1)`, generating
`n[1] × … × n[D]` points in the unit hypercube `[0, 1]^D`.

The points are returned as an array, of size `n`, containing `SVector` of dimension `D`.
"""
Memoize.@memoize Dict function uniform_points(n)
    map(CartesianIndices(n)) do I
        Itup = Tuple(I)
        return SVector((Itup .- 1) ./ (n .- 1))
    end
end
function uniform_points(n, lc, hc)
    @assert length(n) == length(lc) == length(hc) "Dimensions of n, lc, and hc must match."
    return map(x -> lc .+ x .* (hc .- lc), uniform_points(n))
end

Memoize.@memoize function chebyshev_points(nt)
    nodes1d = map(n -> [(-cos(k * π / (n - 1)) + 1) / 2 for k in 0:(n-1)], nt)
    it = Iterators.product(nodes1d...)
    return SVector.(it)
end
function chebyshev_points(nt, lb, ub)
    @assert length(nt) == length(lb) == length(ub) "Dimensions of nt, lb, and ub must match."
    pts = chebyshev_points(nt)
    return map(x -> lb .+ x .* (ub .- lb), pts)
end

"""
    berninterp([T=Float64,] vals::Array, lb, ub)

Construct a Bernstein polynomial of that interpolates the values `vals` at the points given
by `uniform_points(size(vals), lb, ub)`, where `lb` and `ub` are the lower and upper
corners of the hyperrectangle on which the polynomial is defined.

The optional type parameter `T` specifies the precision used for computing the polynomial
coefficients.

# Examples

```jldoctest
using StaticArrays
f = (x) -> (1 - x[1])^2 + x[1]^4 + x[2]^5 * x[1]^3
lb = SVector(0.1, -0.3)
ub = SVector(1.2, 1.7)
pts = ImplicitIntegration.uniform_points((5, 6), lb, ub)
vals = f.(pts)
p = ImplicitIntegration.berninterp(vals, lb, ub)
x = lb .+ (ub - lb) .* rand(SVector{2})
f(x) ≈ p(x)

# output

true
```
"""
function berninterp(T::Type, vals::Array, lb, ub)
    @assert ndims(vals) == length(lb) == length(ub) "Dimensions of vals, lb, and ub must match."
    # create reference interpolation matrix if needed
    n = size(vals)
    V = reference_vandermonde_matrix(T, n)
    c = reshape(V \ vec(vals), n)
    p = BernsteinPolynomial(c, lb, ub)
    return p
end
berninterp(vals::Array, lb, ub) = berninterp(Float64, vals, lb, ub)

"""
    berninterp([T=Float64,] f, n, lb, ub)

Construct a Bernstein polynomial of degree `n` that interpolates the function `f` at the
points given by `uniform_points(n, lb, ub)`, where `lb` and `ub` are the lower and upper
corners of the hyperrectangle on which the polynomial is defined. The optional type parameter
`T` specifies the precision used for computing the polynomial coefficients.
"""
function berninterp(T::Type, f, n, lb, ub)
    @assert length(n) == length(lb) == length(ub) "Dimensions of n, lb, and ub must match."
    pts = uniform_points(n, lb, ub)
    vals = f.(pts)
    return berninterp(T, vals, lb, ub)
end
berninterp(f, n, lb, ub) = berninterp(Float64, f, n, lb, ub)

"""
    reference_vandermonde_matrix([T=Float64,] n)

Vandermond matrix for Bernstein interpolation on `uniform_points(n, lb, ub)`, with `lb = (0, … , 0)` and `ub = (1, … , 1)`. The ordering of the basis elements is the same as in the
evaluation of the Bernstein polynomial. The optional type parameter `T` specifies the
element type of the matrix.
"""
Memoize.@memoize Dict function reference_vandermonde_matrix(T, n)
    N = length(n)
    c = zeros(T, n)
    l = ntuple(i -> zero(T), N) |> SVector
    u = ntuple(i -> one(T), N) |> SVector
    pts = uniform_points(n, l, u)
    p = BernsteinPolynomial(c, l, u)
    A = Matrix{T}(undef, length(c), length(c))
    for i in LinearIndices(c)
        c[i] = one(T)
        for j in LinearIndices(pts)
            A[j, i] = p(pts[j])
        end
        c[i] = zero(T)
    end
    return factorize(A)
end
reference_vandermonde_matrix(n) = reference_vandermonde_matrix(Float64, n)

function Base.show(io::IO, ::MIME"text/plain", p::BernsteinPolynomial{N}) where {N}
    lb = low_corner(p)
    ub = high_corner(p)
    print(
        io,
        "Bernestein polynomial of degree ",
        degree(p),
        " on ",
        '[',
        lb[1],
        ',',
        ub[1],
        ']',
    )
    for i in 2:N
        print(io, " × [", lb[i], ',', ub[i], ']')
    end
end

# Polynomial type for convenience

"""
    struct Polynomial{D,T}

`D`-dimensional polynomial with coefficients of type `T`. The coefficients are stored as a
dense array, and are implicitly associated with a monomial basis; i.e. the `c[I]`
coefficient multiplies the monomial term `prod(x.^(I .- 1))`, where `I` is a `D`-dimensional
multi-index.

Passing a `Dict{NTuple{N,Int},T}` to the constructor will create a polynomial with
coefficients given by `c[k]` for the monomial `x₁^k[1] * x₂^k[2] * ... * x_N^k[N]`, where
`k` is a multi-index of length `N` corresponding to the key, and `c::T` is the value
associated with that key.
"""
struct Polynomial{D,T}
    coeffs::Array{T,D}
end

coefficients(p::Polynomial) = p.coeffs

function Polynomial(c::Dict{NTuple{N,Int},T}) where {N,T}
    degree = ntuple(i -> 0, N)
    for k in keys(c)
        degree = max.(k, degree)
    end
    c′ = zeros(T, degree .+ 1)
    for (k, v) in c
        c′[(k .+ 1)...] = v
    end
    return Polynomial(c′)
end

function Base.show(io::IO, ::MIME"text/plain", p::Polynomial{D}) where {D}
    str = ""
    C = coefficients(p)
    for I in CartesianIndices(C)
        c = C[I]
        iszero(c) && continue
        θ = Tuple(I) .- 1
        x_str = prod(d -> if iszero(θ[d])
            ""
        else
            "x$d^$(θ[d])"
        end, 1:D)
        times_str = sum(θ) == 0 ? "" : " * "
        str = str * "$c $times_str " * x_str * " + "
    end
    return print(io, str[1:(end-2)])
end

"""
    BernsteinPolynomial(p::Polynomial, lb, ub)

Convert a `Polynomial` `p` with coefficients in the monomial basis to a Bernstein polynomial
on the domain defined by the lower corner `lb` and upper corner `ub`.
"""
function BernsteinPolynomial(p::Polynomial, lb, ub)
    c = coefficients(p)
    b = _monomial_to_bernstein(c, lb, ub)
    return BernsteinPolynomial(b, lb, ub)
end

"""
    _monomial_to_bernstein(c::Array, lb, ub)

Convert a polynomial with coefficients `c` in the monomial basis to a Bernstein polynomial
on the hyperrectangle defined by `lb` and `ub`.
"""
function _monomial_to_bernstein(c, lb, ub)
    D = ndims(c)
    b = zero(c)
    c = _rebase(c, lb, ub)
    k = size(c) .- 1
    for i in CartesianIndices(c)
        temp = zeros(Tuple([i[j] for j in 1:D]))
        for l in CartesianIndices(temp)
            temp[l] =
                c[l] * prod(1:D) do j
                    return binomial(i[j] - 1, l[j] - 1) / binomial(k[j], l[j] - 1)
                end
        end
        b[i] = sum(temp)
    end
    return b
end

function _rebase(a::Vector{<:Real}, l::Real, r::Real)
    n = length(a)
    ã = copy(a)
    for i in 0:(n-2)
        ã[(n-i):n] .*= (r - l)
        ã[(n-i-1):(n-1)] .+= ã[(n-i):n] ./ (r - l) .* l
    end
    return ã
end

function _rebase(A::Array{<:Real,D}, L, R) where {(D)}
    Ã = copy(A)
    for d in 1:D
        Ã = mapslices(Ã; dims = d) do a
            return _rebase(a, L[d], R[d])
        end
    end
    return Ã
end

## interface methods for `integrate` and `quadgen`
function bound(b::BernsteinPolynomial, lc, hc)
    if lc == low_corner(b) && hc == high_corner(b)
        return bound(b)
    else
        error("Bounding box does not match polynomial bounds.")
    end
end

function gradient(p::BernsteinPolynomial{N}) where {N}
    return ntuple(d -> derivative(p, d), N)
end

function (∇b::NTuple{N,<:BernsteinPolynomial})(x::SVector{N}) where {N}
    return ntuple(i -> ∇b[i](x), N) |> SVector
end

function bound(∇b::NTuple{N,<:BernsteinPolynomial}, lc, hc) where {N}
    if lc == low_corner(∇b[1]) && hc == high_corner(∇b[1])
        return ntuple(i -> bound(∇b[i]), N) |> SVector
    else
        error("Bounding box does not match polynomial bounds.")
    end
end

function project(p::BernsteinPolynomial{N}, k::Int, v) where {N}
    lc, hc = low_corner(p), high_corner(p)
    if v ≈ lc[k]
        return lower_restrict(p, k)
    elseif v ≈ hc[k]
        return upper_restrict(p, k)
    else
        error("Projection value $v does not match the bounds of the polynomial.")
    end
end

function split(p::BernsteinPolynomial, lb, ub, dir)
    if lb == low_corner(p) && ub == high_corner(p)
        split(p, dir, 0.5)
    else
        error("Bounding box does not match polynomial bounds.")
    end
end

"""
    bernstein_interp(vals, pts, lb, ub)

Construct a Bernstein polynomial on `(lb[1], ub[1]) × … × (lb[N], ub[N])` that interpolates
the values `vals` at the points `pts`. Note that `vals` and `pts`...
"""
function bernstein_interp(vals, pts, lb, ub, degree = size(vals) .- 1)
    return V = vandermonde_matrix(degree, pts, lb, ub)
end

function vandermonde_matrix(degree, pts, lb, ub)
    c = zeros(degree .+ 1)
    p = BernsteinPolynomial(c, lb, ub)
    A = Matrix{Float64}(undef, length(c), length(c))
    for i in LinearIndices(c)
        c[i] = 1.0
        for j in LinearIndices(pts)
            A[j, i] = p(pts[j])
        end
        c[i] = 0.0
    end
    return A
end

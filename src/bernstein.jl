struct BernsteinPolynomial{N,T} <: Function
    coeffs::Array{T,N}
    low_corner::SVector{N,T}
    high_corner::SVector{N,T}
end

function BernsteinPolynomial(c::AbstractArray, lc, hc)
    N = ndims(c)
    return BernsteinPolynomial{N,eltype(c)}(c, SVector{N}(lc), SVector{N}(hc))
end

coefficients(p::BernsteinPolynomial) = p.coeffs
low_corner(p::BernsteinPolynomial)   = p.low_corner
high_corner(p::BernsteinPolynomial)  = p.high_corner
degree(p::BernsteinPolynomial)       = size(coefficients(p)) .- 1

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

# evaluation
function (p::BernsteinPolynomial{N})(x::SVector{N}) where {N}
    l = low_corner(p)
    r = high_corner(p)
    x₀ = (x - l) ./ (r - l)
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

# derivative
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

# split
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
                @. b[1:n-1] = b[1:n-1] * (1 - α) + b[2:n] * α
                b[n] = b[n] * (1 - α)
            else
                push!(c1, b[1])
                pushfirst!(c2, b[i+1])
                @. b[1:i] = b[1:i] * (1 - α) + b[2:i+1] * α
            end
        end
        push!(c1, b[1])
        append!(c1, c2)
        return c1
    end
    lc, hc = low_corner(p), high_corner(p)
    split_pos = lc[d] + (hc[d] - lc[d]) * α
    p1 = BernsteinPolynomial(
        collect(selectdim(coeffs, d, 1:k+1)),
        lc,
        setindex(hc, split_pos, d),
    )
    p2 = BernsteinPolynomial(
        collect(selectdim(coeffs, d, k+1:k+n)),
        setindex(lc, split_pos, d),
        hc,
    )
    return p1, p2
end

# restriction
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

# bound
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

function vandermonde_matrix(degree, pts)
    N = length(degree)
    c = zeros(degree .+ 1)
    l = ntuple(i -> 0.0, N) |> SVector
    u = ntuple(i -> 1.0, N) |> SVector
    p = BernsteinPolynomial(c, l, u)
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
function vandermond_matrix(degree)
    # uniform grid on [0,1]^N for the interpolation
    pts = map(CartesianIndices(degree .+ 1)) do I
        Itup = Tuple(I)
        return SVector((Itup .- 1) ./ degree)
    end
    return vandermonde_matrix(degree, pts)
end

# simple point distributions for interpolation
function uniform_points(degree, lc, hc)
    map(CartesianIndices(degree .+ 1)) do I
        Itup = Tuple(I)
        return SVector((Itup .- 1) ./ degree) .* (hc .- lc) .+ lc
    end
end

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

# add a Polynomial type for convenience

"""
    struct Polynomial{D,T}

`D`-dimensional polynomial with coefficients of type `T`. The coefficients are
stored as a dense array, and are implicitly associated with a monomial basis;
the `c[I]` coefficient multiplies the monomial term `prod(x.^(I .- 1))`, where
`I` is a `D`-dimensional multi-index.
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

function rebase(a::Vector{<:Real}, l::Real, r::Real)
    n = length(a)
    ã = copy(a)
    for i in 0:n-2
        ã[n-i:n] .*= (r - l)
        ã[n-i-1:n-1] .+= ã[n-i:n] ./ (r - l) .* l
    end
    return ã
end

function rebase(A::Array{<:Real,D}, L, R) where {(D)}
    Ã = copy(A)
    for d in 1:D
        Ã = mapslices(Ã; dims = d) do a
            return rebase(a, L[d], R[d])
        end
    end
    return Ã
end

"""
    power2bernstein(a::Array{<:Real,D}, U::HyperRectangle{D}=□(D), k=size(a).-1) where{D}

Convert a polynomial in power series into a Bernstein polynomial on `U` of degree `k`.
"""
function BernsteinPolynomial(p::Polynomial, lb, ub)
    c = coefficients(p)
    b = _monomial_to_bernstein(c, lb, ub)
    return BernsteinPolynomial(b, lb, ub)
end

function _monomial_to_bernstein(c, lb, ub)
    D = ndims(c)
    b = zero(c)
    c = rebase(c, lb, ub)
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

## interface methods
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

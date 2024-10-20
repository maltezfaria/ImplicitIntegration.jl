"""
    struct SubFunction{M,N,T,F}

Represents a function `f̃ : ℝᴹ → ℝ` given by restricting `N` coordinates of
another function `f : ℝᴹ⁺ᴺ → ℝ` to fixed values.

# Fields

  - `f::F`: The underlying function.
  - `dims::SVector{N,Int}`: The dimensions to restrict.
  - `vals::SVector{N,T}`: The values on the restricted dimensions.

The type parameter `M` has to be explicitly provided (e.g. `SubFunction{M}(f, dims, vals)`) since it cannot be inferred from the fields.
"""
struct SubFunction{M,N,T,F}
    f::F # underlying function
    dims::SVector{N,Int} # dimensions to restrict
    vals::SVector{N,T} # values on the restricted dimensions
end

function SubFunction{M}(f::F, dims::SVector{N,Int}, vals::SVector{N,T}) where {M,N,T,F}
    return SubFunction{M,N,T,F}(f, dims, vals)
end

function (f̃::SubFunction{M,N,T})(x̃::SVector{M,T}) where {M,N,T}
    x = multi_insert(x̃, f̃.dims, f̃.vals)
    return f̃.f(x)
end

function gradient(f̃::SubFunction, x̃)
    x = multi_insert(x̃, f̃.dims, f̃.vals)
    ∇f = gradient(f̃.f, x)
    return multi_deleteat(∇f, f̃.dims)
end

function bound(f̃::SubFunction, Ũ::HyperRectangle)
    ã, b̃ = bounds(Ũ)
    a, b = multi_insert(ã, f̃.dims, f̃.vals), multi_insert(b̃, f̃.dims, f̃.vals)
    U = HyperRectangle(a, b)
    return bound(f̃.f, U)
end

function bound_gradient(f̃::SubFunction, Ũ::HyperRectangle)
    ã, b̃ = bounds(Ũ)
    a, b = multi_insert(ã, f̃.dims, f̃.vals), multi_insert(b̃, f̃.dims, f̃.vals)
    U = HyperRectangle(a, b)
    return multi_deleteat(bound_gradient(f̃.f, U), f̃.dims)
end

function restrict(f::SubFunction{M}, U::HyperRectangle{M}, k) where {M}
    lc, hc = bounds(U)
    # idx = findfirst(i -> k > i, f.dims)
    idx = length(f.dims) + 1
    dims = insert(f.dims, idx, k)
    vals_lower = insert(f.vals, idx, lc[k])
    vals_upper = insert(f.vals, idx, hc[k])
    fl = SubFunction{M - 1}(f.f, dims, vals_lower)
    fu = SubFunction{M - 1}(f.f, dims, vals_upper)
    return fl, fu
end

function restriction_type(::Type{SubFunction{M,N,T,F}}) where {M,N,T,F}
    return SubFunction{M - 1,N + 1,T,F}
end

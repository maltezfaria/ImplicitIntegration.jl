"""
    struct HyperRectangle{N,T<:AbstractFloat}

A struct representing a hyperrectangle in N-dimensional space.

# Fields
- `lc::SVector{N,T}`: The lower corner of the hyperrectangle.
- `hc::SVector{N,T}`: The upper corner of the hyperrectangle.

"""
struct HyperRectangle{N,T<:AbstractFloat}
    lc::SVector{N,T}
    hc::SVector{N,T}
end

"""
    const Segment{T} = HyperRectangle{1,T}

A one-dimensional hyperrectangle.
"""
const Segment{T} = HyperRectangle{1,T}

"""
    bounds(rect::HyperRectangle)

Get the lower and upper bounds of a `HyperRectangle`.

# Arguments
- `rect::HyperRectangle`: The `HyperRectangle` object.

# Returns
- A tuple `(lc, hc)` representing the lower and upper bounds of the `HyperRectangle`.
"""
bounds(rect::HyperRectangle) = (rect.lc, rect.hc)

"""
    remove_dimension(rect::HyperRectangle, k)

Remove a dimension from a `HyperRectangle` by deleting the `k`-th element from the lower and upper corners.

# Arguments
- `rect::HyperRectangle`: The input hyperrectangle.
- `k`: The index of the dimension to be removed.

# Returns
A new `HyperRectangle` with the `k`-th dimension removed.
"""
function remove_dimension(rect::HyperRectangle, k)
    lc, hc = bounds(rect)
    return HyperRectangle(deleteat(lc, k), deleteat(hc, k))
end

"""
    split(U::HyperRectangle, dir)

Split a hyperrectangle `U` along the specified direction `dir`.

# Arguments
- `U::HyperRectangle`: The hyperrectangle to be split.
- `dir`: The direction along which to split the hyperrectangle.

# Returns
- `Uₗ`: The left half of the split hyperrectangle.
- `Uᵣ`: The right half of the split hyperrectangle.
"""
function split(U::HyperRectangle, dir)
    lc, hc = bounds(U)
    mid = (lc[dir] + hc[dir]) / 2
    Uₗ = HyperRectangle(lc, setindex(hc, mid, dir))
    Uᵣ = HyperRectangle(setindex(lc, mid, dir), hc)
    return Uₗ, Uᵣ
end

"""
    multi_insert(x::SVector{N,T}, idxs::SVector{M,Int}, vals::SVector{M,T}, ::Val{CHECK})

Like `StaticArrays.insert`, but inserts multiple values into `x`.

Passing a sorted `idxs` together with `Val(false)` for `CHECK` can be used to
assert that `idxs` is sorted without actually checking it.
"""
@generated function multi_insert(
    x::SVector,
    idxs::SVector{M,Int},
    vals::SVector{M},
    ::Val{CHECK} = Val(true),
) where {M,CHECK}
    # Create a sequence of nested calls to `insert` without intermediate
    # assignments
    expr = :x
    for i in 1:M
        expr = :(insert($expr, idxs[$i], vals[$i]))
    end
    return quote
        if CHECK && !issorted(idxs)
            p = sortperm(idxs)
            idxs, vals = idxs[p], vals[p]
        end
        $expr
    end
end

@generated function multi_deleteat(
    x::SVector,
    idxs::SVector{M,Int},
    ::Val{CHECK} = Val(true),
) where {M,CHECK}
    expr = :x
    for i in 1:M
        expr = :(deleteat($expr, idxs[$i]))
    end
    return quote
        if CHECK && !issorted(idxs; rev = true)
            idxs = sort(idxs; rev = true)
        end
        $expr
    end
end

@inline svector(f, N) = SVector(ntuple(f, N))

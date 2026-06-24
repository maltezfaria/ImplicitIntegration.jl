"""
    integrate_threaded(f, ϕ, lc, hc; partition = 2, kwargs...) -> (; val, logger)

Multithreaded version of [`integrate`](@ref). The bounding box `[lc, hc]` is split into a
regular grid of subboxes (`partition` subdivisions per dimension, or an `NTuple` of
per-dimension subdivisions), each subbox is integrated independently on a separate task, and
the results are summed.

This is correct because the implicit domain restricted to the box is the disjoint union of
its restrictions to the subboxes (the shared faces have measure zero, for both volume and
surface integrals). The summation is performed in a fixed (grid) order, so the returned value
is *deterministic and independent of the number of threads*.

`partition` controls the number of tasks (`prod(partition)`); choose it large enough relative
to `Threads.nthreads()` to balance load (cut cells are more expensive than full/empty ones),
but large enough per subbox to amortize task-spawn overhead.

`integrate` is reentrant on the default interface, so no locking is required. All keyword
arguments of [`integrate`](@ref) are forwarded, except `loginfo`, which is not supported here
(per-task loggers cannot be merged meaningfully) and is ignored with a warning if set. The
requested absolute tolerance `tol` is split over the subboxes (each is integrated to
`tol / prod(partition)`) so that the summed absolute error stays within `tol`.

!!! note

    Splitting can place a subbox boundary at or near a critical point of `ϕ` (a point where
    `∇ϕ = 0`, e.g. the centre of a sphere). The algorithm has no good height direction there
    and falls back to a low-order rule on the resulting small box, so the result may be
    slightly less accurate than the (unsplit) serial `integrate` near such points. The error
    remains within the method's tolerance, but choose `partition` so that boundaries avoid
    known critical points when high accuracy near them is required.

# Examples

```jldoctest; output = false
a, b = (0.0, 0.0), (1.5, 1.5)
ϕ = (x) -> x[1]^2 + x[2]^2 - 1
res = integrate_threaded(x -> 1.0, ϕ, a, b; partition = 4)
res.val ≈ π / 4

# output

true

```
"""
function integrate_threaded(
    f,
    ϕ,
    lc::SVector{N,T},
    hc::SVector{N,T};
    partition::Union{Integer,NTuple{N,<:Integer}} = 2,
    kwargs...,
) where {N,T}
    if get(kwargs, :loginfo, false)
        @warn "`loginfo` is not supported by `integrate_threaded`; ignoring it."
    end
    # forward every kwarg except `loginfo` and `tol` (the latter is rescaled per subbox)
    kw = Base.structdiff(NamedTuple(kwargs), NamedTuple{(:loginfo, :tol)})
    nsplit = partition isa Integer ? ntuple(_ -> Int(partition), N) : map(Int, partition)
    all(>(0), nsplit) || throw(ArgumentError("`partition` must be positive in every dimension"))
    edges = ntuple(d -> range(lc[d], hc[d]; length = nsplit[d] + 1), N)
    cell_ids = collect(Iterators.product(ntuple(d -> 1:nsplit[d], N)...))
    # Split the absolute tolerance over the subboxes so that the summed absolute error stays
    # within the requested `tol` (each subbox contributes at most `tol / ncells`).
    tol_cell = get(kwargs, :tol, 1e-8) / length(cell_ids)
    RET = typeof(f(lc) * one(T) + f(hc) * one(T))
    results = Vector{RET}(undef, length(cell_ids))
    @sync for (i, c) in enumerate(cell_ids)
        Threads.@spawn begin
            a = SVector(ntuple(d -> T(edges[d][c[d]]), N))
            b = SVector(ntuple(d -> T(edges[d][c[d]+1]), N))
            results[i] = integrate(f, ϕ, a, b; tol = tol_cell, kw...).val
        end
    end
    val = sum(results)
    return (; val, logger = nothing)
end

function integrate_threaded(f, ϕ, lc, hc; kwargs...)
    @assert length(lc) == length(hc) "Lower and upper corners must have the same length."
    N = length(lc)
    T = promote_type(float(eltype(lc)), float(eltype(hc)))
    return integrate_threaded(f, ϕ, SVector{N,T}(lc), SVector{N,T}(hc); kwargs...)
end

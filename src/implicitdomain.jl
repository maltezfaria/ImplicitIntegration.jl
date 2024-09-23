"""
    struct ImplicitDomain

Structure used to represent a domain `V` as the intersection of several
implicitly defined domains `Vᵢ`. The domains `Vᵢ` are defined using a box `U`
and a set of `functions` and `flags` as follows:

- `Vᵢ = {x ∈ U : ϕᵢ(x) > 0}` if `sᵢ = 1`.
- `Vᵢ = {x ∈ U : ϕᵢ(x) < 0}` if `sᵢ = -1`.
- `Vᵢ = {x ∈ U : ϕᵢ(x) ≠ 0}` if `sᵢ = 0`.

where `ϕᵢ = functions[i]` is a function mapping `U` to `ℝ`, and `sᵢ = flags[i]`
is a flag used to define the domain `Vᵢ`.

# Fields
- `box::HyperRectangle{N,T}`: bounding box of the domain `V`.
- `functions::Vector{Function}`: vector of functions defining the domains `Vᵢ`.
- `flags::Vector{Int}`: vector of flags `sᵢ` associated with each function.

"""
struct ImplicitDomain{N,T}
    box::HyperRectangle{N,T}
    functions::Vector
    flags::Vector{Int}
end

```@meta
CurrentModule = ImplicitIntegration
```

# ImplicitIntegration

## [Overview](@id overview)

This package provides an [`integrate`](@ref) function to approximate *volume* and *surface*
integrals over implicitly defined domains in arbitrary dimensions. More specifically, it
allows for the computation of **volume integrals** of the form

```math
    \int_{\phi(\boldsymbol{x}) < 0 \ \cap \ U}  f(\boldsymbol{x}) \, \mathrm{dV},
```

and **surface integrals** of the form

```math
    \int_{\phi(\boldsymbol{x}) = 0 \ \cap \ U}  f(\boldsymbol{x}) \, \mathrm{dS},
```

where ``\phi : \mathbb{R}^d \to \mathbb{R}`` is a level-set function implicitly defining the
surface, and ``U = [a_1, b_1] \times \ldots \times [a_d, b_d]`` is a bounding
[`HyperRectangle`](@ref).

!!! tip "Algorithm"
    The algorithm implemented is based on [saye2015high](@cite), and relies on the ability
    to bound functions and their derivatives over hyperrectangles. Reading the paper is a
    good idea if you want to understand the details and limitations of the method.

## Installation

Since the package is not yet registered, you can install it by running

```julia
using Pkg; Pkg.add("https://github.com/maltezfaria/ImplicitIntegration.jl");
```

## Basic usage

### `integrate`

The main function provided by this package is [`integrate`](@ref), which computes implicit
integrals using an adaptive quadrature method. Here is how you can use it to compute the
volume of a sphere:

```@example overview-example
using ImplicitIntegration, StaticArrays
ϕ = (x) -> x[1]^2 + x[2]^2 + x[3]^2 - 1
f = (x) -> 1
lc = (-1.1, -1.1, -1.1)
hc = (1.1, 1.1, 1.1)
int_volume  = integrate(f, ϕ, lc, hc)
nothing # hide
```

The `integrate` function returns a `NamedTuple` object containing the value of the integral
in `val`, as well as a [`logger`](@ref ImplicitIntegration.LogInfo) object containing
information about the computation:

```@example overview-example
println("Computed volume:  $(int_volume.val)")
println("Error of volume:  $(int_volume.val - 4π/3)")
println(int_volume.logger)
```

To compute a surface integral isntead, simply set the `surface` keyword argument to `true`:

```@example overview-example
int_surface = integrate(f, ϕ, lc, hc; surface = true)
println("Computed surface: $(int_surface.val)")
println("Error of surface: $(int_surface.val - 4π)")
println(int_surface.logger)
```

### `quadgen`

For situations where you need to compute the integrals of multiple functions over *the same
domain*, you can use [`quadgen`](@ref) to generate a quadrature instead. It works similarly
to [`integrate`](@ref), but returns a [`Quadrature`](@ref) object instead of the integral
value. Here is how to create a surface quadrature for a [Cassini
oval](https://en.wikipedia.org/wiki/Cassini_oval) in 3D:

```@example overview-example
using ImplicitIntegration, StaticArrays
using GLMakie, LinearAlgebra
p1 = SVector(-1.0, 0, 0)
p2 = SVector(1.0, 0, 0)
b = 1.1
ϕ = x -> (x - p1) ⋅ (x - p1) * (x - p2) ⋅ (x - p2) - b^2
out = quadgen(ϕ, (-2, -2, -2), (2, 2, 2); order = 5, surface = true)
xx = yy = zz = range(-1.5, 1.5, length = 200)
vv = ϕ.(SVector.(Iterators.product(xx, yy, zz)))
volume(xx, yy, zz, vv, algorithm = :iso, transparency = true, alpha = 0.4, isovalue = 0)
scatter!(out.quad.coords, markersize = 2, color = :red)
current_figure()
```

See the [`quadgen`](@ref) docstrings for more information on the available options.

## Interface methods

As alluded to in the [overview section](@ref overview), the underlying algorithm relies on bounding
the level-set function ``\phi`` as well as its partial derivatives ``\partial_{x_i} \phi``
in order to find an appropriate height function. By default, `ImplicitIntegration.jl` uses a
combination of
[IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) and
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) to perform these task. 

There are at least two situations, however, where it may be beneficial (or even required) to
overload the defaults:

- `ForwardDiff` and/or `IntervalArithmetic` do not work with your function type
- The bounds provided by the default methods are innacurate and/or slow

In such cases, you can overload the following methods for your levelset `ϕ`:

- `ϕ(x::SVector{N,<:Real})`: evaluate the function
- `ϕ(I::SVector{N,<:Interval})`: bound the function
- `ϕ(x::SVector{N,<:ForwardDiff.Dual{Tg,T<:Real,M}})`, `M <= N`: evaluate the Jacobian
  vector product
- `ϕ(x::SVector{N,<:ForwardDiff.Dual{Tg,T<:Interval,M}})`, `M <= N`: bound the Jacobian
  vector product
  
For instance, here's how to implement the Jacobian vector product overload using a
user-provided `ϕ_and_∇ϕ(x)` function:
```
function ϕ(x::SVector{S,<:ForwardDiff.Dual{Tg,T,Npartials}}) where {S, Tg, T, Npartials}
    ϕx, ∇ϕx = ϕ_and_∇ϕ(ForwardDiff.value.(x))
    partials = ntuple(ipart -> sum(∇ϕx .* ForwardDiff.partials.(x,ipart)), Npartials)
    ForwardDiff.Dual{Tg}(ϕx, partials...)
end
```
This is enough, provided `ϕ_and_∇ϕ` supports interval inputs. If it doesn't, just add
a `ϕ_and_∇ϕ(x::SVector{<:Interval})` method that computes the bound and returns it as
an `Interval`. Note that ForwardDiff and IntervalArithmetic are both composable: if
`ϕ` is the composition of two functions `ϕ1` and `ϕ2`, you can use custom overloads 
for either or both of the two functions.

## Bibliography

```@bibliography
```

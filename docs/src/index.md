```@meta
CurrentModule = ImplicitIntegration
```

# ImplicitIntegration

## Overview

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

## Bibliography

```@bibliography
```

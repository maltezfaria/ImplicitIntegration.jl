# ImplicitIntegration

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://maltezfaria.github.io/ImplicitIntegration.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://maltezfaria.github.io/ImplicitIntegration.jl/dev/)
[![Build Status](https://github.com/maltezfaria/ImplicitIntegration.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/maltezfaria/ImplicitIntegration.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/maltezfaria/ImplicitIntegration.jl/graph/badge.svg?token=0bVypWcVTc)](https://codecov.io/gh/maltezfaria/ImplicitIntegration.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

**ImplicitIntegration** is a *Julia* library for approximating volume and surface integrals
over implicitly defined domains in arbitrary dimensions. It implements the algorithm
described in

> Saye, R. I. (2015). High-order quadrature methods for implicitly defined surfaces and
> volumes in hyperrectangles. SIAM Journal on Scientific Computing, 37(2), A993-A1019.

which allows for the computation of **volume integrals** of the form

```math
    \int_{\phi(\boldsymbol{x}) < 0 \ \cap \ U}  f(\boldsymbol{x}) \, \mathrm{dV},
```

and **surface integrals** of the form

```math
    \int_{\phi(\boldsymbol{x}) = 0 \ \cap \ U}  f(\boldsymbol{x}) \, \mathrm{dS},
```

where $\phi : \mathbb{R}^d \to \mathbb{R}$ is a level-set function implicitly defining the
surface, and $U = [a_1, b_1] \times \ldots \times [a_d, b_d]$ is a bounding
box.

Key features include:

- Implicit volume and surface integrals in arbitrary dimensions
- Flexible domain definition through native Julia functions
- High-precision integration
- Customizable interface

## Quick Start

To install, open the Julia REPL and run:

```julia
using Pkg
Pkg.add("ImplicitIntegration")
```

Here is how to compute the area of a unit disk in 2D:

```julia
using ImplicitIntegration
ϕ = (x) -> sqrt(x[1]^2 + x[2]^2) - 1 # level-set function
f = (x) -> 1                         # integrand
result = integrate(f, ϕ, (-1.5, -1.5), (1.5, 1.5))
result.val ≈ π # true
```

To compute the perimeter, simply set `surface = true`:

```julia
result = integrate(f, ϕ, (-1.5, -1.5), (1.5, 1.5); surface = true)
result.val ≈ 2π # true
```

More interesting surfaces can be obtained by changing `ϕ`. Here is an example of a
quadrature for a three-dimensional [Cassini
oval](https://en.wikipedia.org/wiki/Cassini_oval):

```julia
using ImplicitIntegration, StaticArrays, GLMakie, LinearAlgebra
p1, p2, p3 = SVector(-1.0, 0, -√3 / 2), SVector(1.0, 0, -√3 / 2), SVector(0, 0, √3 / 2)
b = 1.5
# surface s.t. product of distances to p1, p2, p3 is constant
ϕ = x -> (x - p1) ⋅ (x - p1) * (x - p2) ⋅ (x - p2) * (x - p3) ⋅ (x - p3) - b^2 
out = quadgen(ϕ, (-2, -2, -2), (2, 2, 2); order = 1, surface = true)
scatter(out.quad.coords; markersize = 5, color = :red)
```

![trifoil](trifoil.png)

For detailed usage, examples, and advanced options, please refer to the [official
documentation](https://maltezfaria.github.io/ImplicitIntegration.jl/dev/).

## Related packages

- [Algoim](https://algoim.github.io): C++ library by the author of the original paper, Robert I. Saye
- [Algoim.jl](https://github.com/ericneiva/Algoim.jl): a Julia wrapper around the C++
  library above
- [QuadratureOnImplicitRegions.jl](https://github.com/HMegh/QuadratureOnImplicitRegions.jl):
  another implementation of Saye's algorithm in Julia

```@meta
CurrentModule = ImplicitIntegration
```

# ImplicitIntegration

## Overview

This package provides an [`integrate`](@ref) function to approximate *volume*
and *surface* integrals over implicitly defined domains in arbitrary dimensions.
The algorithm is based on [saye2015high](@cite).

More specifically, it allows for the computation of *volume integrals* of the form

```math
    \int_{\phi(\boldsymbol{x}) < 0 \ \cap \ U}  f(\boldsymbol{x}) \, \mathrm{dV},
```

and *surface integrals* of the form

```math
    \int_{\phi(\boldsymbol{x}) = 0 \ \cap \ U}  f(\boldsymbol{x}) \, \mathrm{dS},
```

where ``\phi : \mathbb{R}^d \to \mathbb{R}`` is a level-set function, and ``U =
[a_1, b_1] \times \ldots \times [a_d, b_d]`` is a bounding
[`HyperRectangle`](@ref).

## Installation

Since the package is not yet registered, you can install it by running

```julia
using Pkg; Pkg.add("https://github.com/maltezfaria/ImplicitIntegration.jl");
```

## Basic usage

The main function provided by this package is `integrate`:

```@docs
integrate
```

## Bibliography

```@bibliography
Pages = ["index.md"]
Canonical = false
```

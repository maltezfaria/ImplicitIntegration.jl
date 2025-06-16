```@meta
CurrentModule = ImplicitIntegration
```

# [Specializations](@id specializations)

!!! details "Importance of tight bounds"
    The algorithm for implicit integration is described in [saye2015high](@cite) relies on
    the ability to compute bounds on the implicit function and its derivatives over a
    hyper-rectangle. Bounds on the implicit function are used to determine whether a box is
    empty or full, in which case a subdivision is not necessary. Bounds on the gradient are
    used to determine if a height direction exists for the dimensional reduction step, in
    which case the problem is recursively reduced to a lower dimension. If the algorithm is
    unable to prove that a box is empty, full, or that a height direction exists, it
    subdivides the box and repeats the process on the sub-boxes. The tightness of the bounds
    therefore directly affects the efficiency of the algorithm, as tighter bounds lead to
    fewer subdivisions.

`ImplicitIntegration` provides a generic method for bounding functions and their derivatives
using automatic differentiation via `ForwardDiff` and range estimation with
`IntervalArithmetic`. While this approach has the advantage of being rather generic, the
resulting bounds are often conservative and may lead to unnecessary subdivisions during the
quadrature generation process.

For specific classes of functions, we can leverage additional mathematical properties to
compute tighter bounds more efficiently. These specialized bounds can significantly improve
the performance of implicit integration, reducing both the computational cost and the number
of quadrature points needed for accurate results.

This section describes the interface for implementing custom bounds for specific function
types, and discuss one such specialization already implemented in the package: Bernstein
polynomials.

## [Interface](@id interface)

To implement specialized bounds for a specific function type, you need to overload the
following methods:

- [`bound(f::MyType, lc, hc) --> (fmin, fmax)`](@ref ImplicitIntegration.bound): Compute
  bounds for the function `f` over the hyperrectangle defined by lower corner `lc` and upper
  corner `hc`.
- [`gradient(f::MyType) --> ∇f`](@ref ImplicitIntegration.gradient): Compute the gradient of
  the function `f` as a function. The returned object `∇f` should support evaluation on the
  form `∇f(x)` where `x` is a point in the domain of `f`, and `bound(∇f, lc, hc)`.
- [`project(f, k, v)`](@ref ImplicitIntegration.project): Project the function `f` onto the
  `k`-th coordinate using the value `v`. The returned object should support `bound`,
  `gradient`, and `project` methods, allowing it to be used in the same way as the original
  function.
- [`split(f, lb, ub, dir) --> fₗ, fᵤ`](@ref ImplicitIntegration.split): Split the function
  `f` along the specified direction `dir` into two parts, one for the lower half and one for
  the upper half of the hyperrectangle defined by `lb` and `ub`. The returned objects should
  also support `bound`, `gradient`, and `project` methods.

!!! tip "Disabling default implementations"
    When implementing the interface for a specific function type, you may want to disable
    the default implementations provided by `ImplicitIntegration` to make sure you have
    correctly overloaded the proper methods. You can do this by calling
    [`disable_default_interface()`](@ref ImplicitIntegration.disable_default_interface). To
    re-enable the default interface, use [`enable_default_interface()`](@ref
    ImplicitIntegration.enable_default_interface).

## [Bernstein polynomials](@id bernstein-polynomials)

Because the Bernstein basis forms a partition of unity, computing the upper and lower bounds
of a Bernstein polynomial is particularly simple and only requires evaluating the `extrema`
of the Bernstein coefficients. Such a bound is usually much tighter than the bounds obtained
by interval arithmetic. `ImplicitIntegration` provides a `BernsteinPolynomial` type that
implements the interface described above; see the `bernstein.jl` file in `src` for the
implementation details.

There are several ways to construct a `BernsteinPolynomial` object. The most direct one is
to directly provide an array containing the Bernstein coefficients:

```@example bernstein
using ImplicitIntegration
lb = (-1.0, -1.0)
ub = (1.0, 1.0)
c  = rand(3,3) # Bernstein coefficients
p = ImplicitIntegration.BernsteinPolynomial(c, lb, ub)
```

The polynomial can now be evaluated at a point `x`

```@example bernstein
x = (0.5, 0.5)
p(x)
```

or bounded over its domain:

```@example bernstein
ImplicitIntegration.bound(p)
```

Another way to create a `BernsteinPolynomial` is to interpolate a function over the grid
points using [`berninterp`](@ref), as illustrated below:

```@example bernstein
lb = (-1.0, -1.0)
ub = (1.0, 1.0)
npts = (4,4)
interp_pts = ImplicitIntegration.uniform_points(npts, lb, ub)
f = x -> x[1]^2 + x[2]^2 - 1.0 # implicit function
p = ImplicitIntegration.berninterp(f.(interp_pts), lb, ub)
```

Since in this example `f` is itself a polynomial, and we interpolate it using a polynomial
of higher degree, the interpolant should be exact (up to rounding errors):

```@example bernstein
using GLMakie
@assert all(abs(p((x,y)) - f((x,y))) < 1e-14 for x in lb[1]:0.1:ub[1], y in lb[2]:0.1:ub[2]) # hide
xx = lb[1]:0.05:ub[1]
yy = lb[2]:0.05:ub[2]
vals = [p((x,y)) - f((x,y)) for x in lb[1]:0.1:ub[1], y in lb[2]:0.1:ub[2]]
fig, ax, hm = heatmap(xx, yy, vals)
Colorbar(fig[1,2], hm; label = "|p(x) - f(x)|", labelrotation = 0)
current_figure() # hide
```

!!! warning "Increasing the interpolation degree"
    The `berninterp` function illustrated above performs polynomial interpolating using a
    Bernstein basis on a uniform grid. As is well known, the Lebesgue constant for such
    interpolants grows exponentially with the polynomial degree. In this package you usually
    want to use the `berninterp` function with a fixed degree, and decrease *mesh size*
    parameter to obtain convergence. An alternative, not yet implemented, is to interpolate
    the function on a Chebyshev grid, and then convert the resulting Chebyshev polynomial
    into the Bernstein form.

Finally, `ImplicitIntegration` can automatically handle polynomials constructed using the
`DynamicPolynomials.jl` package. When you use types from this package, the implicit function
is automatically converted to a `BernsteinPolynomial` for bounds computation. This allows
for efficient integration of polynomial implicit functions without requiring manual
conversion to the Bernstein form. This specialization provides tighter bounds for the
function value and its gradient, helping the recursive subdivision algorithm to stop
earlier. The following example illustrates this difference in practice:

```@example polynomials
using ImplicitIntegration, StaticArrays
using GLMakie, LinearAlgebra, DynamicPolynomials
p1 = SVector(-1.0, 0)
p2 = SVector(1.0, 0)
b = 1.01
ϕ = x -> (x - p1) ⋅ (x - p1) * (x - p2) ⋅ (x - p2) - b^2
@polyvar x y
poly = ((x,y) .- p1) ⋅ ((x,y) .- p1) * ((x,y) .- p2) ⋅ ((x,y) .- p2) - b^2
quad, logger = quadgen(ϕ, (-2, -2), (2, 2); order = 1, surface = true, loginfo=true)
quad_poly, logger_poly = quadgen(poly, (-2, -2), (2, 2); order = 1, surface = true, loginfo=true)
## plot
fig = Figure()
ax1 = Axis(fig[1, 1], aspect = DataAspect(), title = "Default bounds", xlabel = "x", ylabel = "y")
scatter!(ax1, quad.coords, markersize = 10, color = :red)
plot!(ax1, logger)
ax2 = Axis(fig[1, 2], aspect = DataAspect(), title = "Bernstein bounds", xlabel = "x", ylabel = "y")
scatter!(ax2, quad_poly.coords, markersize = 10, color = :red)
plot!(ax2, logger_poly)
current_figure() # hide
```

The figure above illustrates the difference between using default bounds (based on interval
arithmetic) and Bernstein polynomial bounds for quadrature generation (triggered by passing
a `DynamicPolynomial.Poly` object to `quadgen`). The right panel shows how using Bernstein
polynomial bounds results in more efficient and accurate quadrature points (red dots) around
the zero level set of the implicit function. Note the significant reduction in the number of
evaluated boxes and the tighter approximation of the curve.

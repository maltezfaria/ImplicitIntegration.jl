```@meta
CurrentModule = ImplicitIntegration
```

# [Specializations](@id specializations)

## [Interface](@id interface)

## [Bernstein polynomials](@id bernstein-polynomials)

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

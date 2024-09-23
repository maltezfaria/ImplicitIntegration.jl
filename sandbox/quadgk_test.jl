using QuadGK
using CairoMakie

pts = Float64[]
f = x -> begin
    push!(pts, x)
    abs(x)
end

I, E = quadgk(f, -π, sqrt(2); atol = 1e-4)

@info "Integral: $I, error: $E, num of pts: $(length(pts))"

scatter(pts, 0 * pts; markersize = 5, color = :blue)

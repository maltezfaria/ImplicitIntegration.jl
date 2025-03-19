using ImplicitIntegration
using Test
using StaticArrays
using GLMakie

function plot_hyperrectangle!(ax, h::ImplicitIntegration.HyperRectangle)
    lc = h.lc
    hc = h.hc
    x = [lc[1], hc[1], hc[1], lc[1], lc[1]]
    y = [lc[2], lc[2], hc[2], hc[2], lc[2]]
    return lines!(ax, x, y; color = :black)
end

function plot_tree(root::ImplicitIntegration.TreeNode)
    fig = Figure()
    ax = Axis(fig[1, 1]; aspect = DataAspect())
    plot_tree!(ax, ImplicitIntegration.lift(root))
    return fig
end
function plot_tree!(ax, root::ImplicitIntegration.TreeNode)
    isempty(root.children) && plot_hyperrectangle!(ax, root.box)
    for (child, dir) in root.children
        @assert dir == 0
        plot_tree!(ax, child)
    end
    return nothing
end

##
using ImplicitIntegration
using Test
using StaticArrays
using GLMakie
δ = 0
a, b = SVector(0.0, 0.0), SVector(1.0 + δ, 1.0 + δ)
ϕ = (x) -> sqrt(x[1]^2 + x[2]^2) - 1
I, loginfo = ImplicitIntegration.integrate(x -> 1.0, ϕ, a, b; loginfo = true)
plot(loginfo)

# plot_tree(loginfo.tree)
# contour!(
#     a[1]:0.01:b[1],
#     a[2]:0.01:b[2],
#     (x, y) -> ϕ(SVector(x, y));
#     levels = [0.0],
#     linewidth = 2,
#     color = :red,
# )
# current_figure()

##
using ImplicitIntegration
using Test
using StaticArrays
using GLMakie
δ = 0.1
a, b = SVector(0.0, 0.0, 0.0), SVector(1.0 + δ, 1.0 + δ, 1.0 + δ)
ϕ = (x) -> sqrt(x[1]^2 + x[2]^2 + x[3]^3) - 1
I, loginfo = ImplicitIntegration.integrate(x -> 1.0, ϕ, a, b; loginfo = true)
plot(loginfo)

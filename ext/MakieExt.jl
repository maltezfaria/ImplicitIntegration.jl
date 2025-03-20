module MakieExt

using Makie
import ImplicitIntegration

function __init__()
    @info "Loading Makie extension for ImplicitIntegration.jl"
end

Makie.@recipe(TreePlot, loginfo) do scene
    return Theme()
end

function Makie.plot!(p::TreePlot)
    loginfo = p.loginfo
    tree = @lift ($loginfo).tree
    N = @lift ndims($tree)
    ImplicitIntegration.lift_all_dims(to_value(tree))
    if to_value(N) == 2
        _plot_tree_2d(p, to_value(tree))
    elseif to_value(N) == 3
        _plot_tree_3d(p, to_value(tree))
    else
        throw(ArgumentError("Plot of $(to_value(N)) dimensional tree not supported."))
    end
    return p
end

# function plot_hyperrectangle!(ax, h::ImplicitIntegration.HyperRectangle)
#     lc = h.lc
#     hc = h.hc
#     x = [lc[1], hc[1], hc[1], lc[1], lc[1]]
#     y = [lc[2], lc[2], hc[2], hc[2], lc[2]]
#     return lines!(ax, x, y; color = :black)
# end

function _plot_tree_2d(p, node)
    if isempty(node.children)
        lc = node.box.lc
        hc = node.box.hc
        x = [lc[1], hc[1], hc[1], lc[1], lc[1]]
        y = [lc[2], lc[2], hc[2], hc[2], lc[2]]
        return lines!(p, x, y; color = :black)
    else
        for (child, dir) in node.children
            @assert dir == 0
            _plot_tree_2d(p, child)
        end
    end
end

function _plot_tree_3d(p, node)
    if isempty(node.children)
        lc = node.box.lc
        hc = node.box.hc
        vertices = [
            lc,
            [hc[1], lc[2], lc[3]],
            [hc[1], hc[2], lc[3]],
            [lc[1], hc[2], lc[3]],
            [lc[1], lc[2], hc[3]],
            [hc[1], lc[2], hc[3]],
            hc,
            [lc[1], hc[2], hc[3]],
        ]
        edges = [
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 1), # Bottom face
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 5), # Top face
            (1, 5),
            (2, 6),
            (3, 7),
            (4, 8),  # Vertical edges
        ]
        for (i, j) in edges
            lines!(
                p,
                [vertices[i][1], vertices[j][1]],
                [vertices[i][2], vertices[j][2]],
                [vertices[i][3], vertices[j][3]];
                color = :black,
            )
        end
    else
        for (child, dir) in node.children
            @assert dir == 0
            _plot_tree_3d(p, child)
        end
    end
end

Makie.plottype(::ImplicitIntegration.LogInfo) = TreePlot

## Pick correct axis type based on dimension
function Makie.preferred_axis_type(p::TreePlot)
    loginfo = p.loginfo
    tree = @lift $loginfo.tree
    dim = @lift ndims($tree)
    if to_value(dim) == 2
        return Axis
    elseif to_value(dim) == 3
        return LScene
    else
        throw(ArgumentError("Plot of $dim dimensional tree not supported."))
    end
end

end # module

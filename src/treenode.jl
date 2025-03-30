"""
    mutable struct TreeNode
"""
mutable struct TreeNode
    box::HyperRectangle
    # each children is either a TreeNode of the same dimension as the parent, in which case
    # `dir` is 0, or a TreeNode of one dimension less, in which case `dir` is the dimension
    # of the parent that remains unchanged.
    children::Vector{Tuple{TreeNode,Int}}
end
Base.ndims(node::TreeNode) = ndims(node.box)

TreeNode(box::HyperRectangle) = TreeNode(box, [])

"""
    lift_all_dims(root::TreeNode)

Lift all nodes in `root` to the same dimension as `root`.
"""
function lift_all_dims(root::TreeNode)
    n = ndims(root)
    for dim in 1:n-1
        lift_dim(root, dim)
    end
    return root
end

"""
    lift_dim(node, dim)

Lift all boxes of dimension `dim` in the tree rooted at `node` to the same dimension as their parent.
"""
function lift_dim(node, dim)
    for i in eachindex(node.children)
        child, dir = node.children[i]
        if ndims(child) == dim
            lift_node_and_children(child, dir, node.box.lc[dir], node.box.hc[dir])
            node.children[i] = (child, 0)
        else
            lift_dim(child, dim)
        end
    end
    return node
end

"""
    lift_node_and_children(node::TreeNode, dim, lb, ub)

Lift each `box` in the tree rooted at `node` by appending `lb` and `ub` to the `dim`-th
dimension of the `box`.
"""
function lift_node_and_children(node::TreeNode, k, lb, ub)
    node.box = HyperRectangle(insert(node.box.lc, k, lb), insert(node.box.hc, k, ub))
    for (child, dir) in node.children
        @assert dir == 0
        lift_node_and_children(child, k, lb, ub)
    end
end

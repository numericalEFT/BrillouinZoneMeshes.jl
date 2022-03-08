using AbstractTrees

struct MyNode
    val::Int
    children::Vector{MyNode}
end

AbstractTrees.children(node::MyNode) = node.children

function MyNode(val::Int)
    val = val
    children = Vector{MyNode}()
    for i in 1:val
        push!(children, MyNode(val-1))
    end

    return MyNode(val, children)
end

tree = MyNode(3)

print_tree(tree)


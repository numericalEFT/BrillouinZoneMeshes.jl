using SpaceGrid
using SpaceGrid.AbstractTrees, SpaceGrid.GridTree, SpaceGrid.BaseMesh

mesh = UniformMesh{2, 4}([0.0, 0.0], [0 1 ; 1 0])
show(mesh)

naiveisfine(depth, pos) = depth >= 2

tree = GridNode{2}(naiveisfine)
print_tree(tree)

for node in PostOrderDFS(tree)
    if isempty(node.children)
        println(node.pos)
    end
end

tg = uniformtreegrid(naiveisfine, [0 1; 1 0])

println(tg.subgrids)

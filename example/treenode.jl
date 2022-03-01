
using SpaceGrid: AbstractTrees, GridTree

naiveisfine(depth, pos) = depth >= 2

AbstractTrees.print_tree(GridTree.GridNode{2}(naiveisfine))


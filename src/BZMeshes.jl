module BZMeshes

using ..StaticArrays
using ..LinearAlgebra

using ..AbstractMeshes
using ..AbstractMeshes: _inds2ind, _ind2inds
using ..Cells
using ..BaseMesh
using ..MeshMaps
using ..PointSymmetry
using ..CompositeGrids
using ..Roots
# using ..BaseMesh.AbstractUniformMesh
import ..showfieldln
import ..showfield
import ..SHOWINDENTION

export UniformBZMesh, DFTK_Monkhorst_Pack, Monkhorst_Pack

struct OnBrillouin <: MeshDomain end # has cell::Cell as 

AbstractMeshes.lattice_vector(::OnBrillouin, mesh) = mesh.cell.recip_lattice
AbstractMeshes.inv_lattice_vector(::OnBrillouin, mesh) = mesh.cell.inv_recip_lattice
AbstractMeshes.lattice_vector(::OnBrillouin, mesh, i::Int) = Cells.get_latvec(mesh.cell.recip_lattice, i)
AbstractMeshes.inv_lattice_vector(::OnBrillouin, mesh, i::Int) = Cells.get_latvec(mesh.cell.inv_recip_lattice, i)
AbstractMeshes.cell_volume(::OnBrillouin, mesh) = mesh.cell.recip_cell_volume

"""
    struct UniformBZMesh{T, DIM} <: AbstractUniformMesh{T, DIM}

Uniformly distributed Brillouin zone mesh. Defined as a uniform mesh on 1st Brillouin zone
with Brillouin zone information stored in mesh.cell::Cell. 

# Parameters:
- `T`: type of data
- `DIM`: dimension of the Brillouin zone

# Members:
- `cell`: Cell information including lattice info, atom pos and allowed G vectors
- `origin`: origin of the uniform mesh, related to convention of 1st Brillouin zone. Commonly set to either (0,0,0) or such that (0,0,0) is at the center
- `size`: size of the uniform mesh. For Monkhorst-Pack mesh require even number.
- `shift`: k-shift of each mesh point. Set all to be zero for Gamma-centered and all to be 1//2 for M-P mesh
"""
# struct UniformBZMesh{T,DIM} <: AbstractMesh{T,DIM}
struct UniformBZMesh{T,DIM} <: AbstractUniformMesh{T,DIM}
    cell::Cell{T,DIM}

    origin::SVector{DIM,T}
    size::NTuple{DIM,Int}
    shift::SVector{DIM,Rational}
end

# default shift is 1/2, result in Monkhorst-Pack mesh
# with shift = 0, result in Gamma-centered
# can also customize with shift::SVector by calling default constructor
# \Gamma=(0,0,0) is at center by default, can be set at corner by setting origin to it
"""
    function UniformBZMesh(; br::Cell, origin, size, shift)

customized constructor for UniformBZMesh. The parameters origin and shift is provided to customize
the mesh as Gamma-centered or M-P mesh. 

# Parameters:
- `cell`: Cell info
- `origin`: a number indicating shift of origin. 
    the actuall origin becomes origin*(b1+b2+b3)
    default value origin=-1/2 takes (0,0,0) to center of 1st BZ, origin=0 makes mesh[1,1,1]=(0,0,0)+shift
- `size`: size of the mesh
- `shift`: additional k-shift for mesh points. 
    actuall shift is shift*(b1/N1+b2/N2+b3/N3)
    for even N, shift=1/2 avoids high symmetry points while preserve symmetry.
"""
function UniformBZMesh(;
    cell::Cell{T,DIM},
    origin::Real=-1 // 2,
    size::Union{AbstractVector,Tuple},
    shift::Union{AbstractVector{Bool},AbstractVector{Int}}=[true for _ in eachindex(size)]
) where {T,DIM}

    _shift = [s == true ? 1 // 2 : 0 for s in shift]

    return UniformBZMesh{T,DIM}(
        cell, (cell.recip_lattice * ones(T, DIM)) * origin, tuple(size...), _shift
    )
end

function DFTK_Monkhorst_Pack(;
    cell::Cell{T,DIM},
    size,
    shift::AbstractVector{Bool}=[false, false, false]) where {T,DIM}
    kshift = [iseven(size[i]) ? shift[i] : !(shift[i]) for i in 1:DIM]
    return UniformBZMesh(
        cell=cell,
        origin=-1 // 2, #origin
        size=tuple(size...),
        shift=kshift
    )
end

function Monkhorst_Pack(;
    # Gamma_centered: origin=0, 
    # Monkhorst-Pack: origin=-1/2, consistent with VASP
    # to be consistent with DFTK: 
    #  - N is even, VASP is the same as DFTK: shift=0 will include Gamma point, shift=1/2 will not
    #  - N is odd, VASP is different as DFTK: shift=0 will not include Gamma point, shift=1/2 will
    cell::Cell{T,DIM},
    size,
    shift::AbstractVector=[false, false, false]
) where {T,DIM}
    # kshift = [(iseven(size[i]) ? shift[i] : shift[i] + 1 // 2) for i in 1:DIM]
    return UniformBZMesh(
        cell=cell,
        origin=-1 // 2,
        size=tuple(size...),
        shift=shift
    )
end

AbstractMeshes.MeshDomain(::Type{<:UniformBZMesh}) = OnBrillouin()

# function Base.show(io::IO, mesh::UniformBZMesh)
#     println("UniformBZMesh with $(length(mesh)) mesh points")
# end

function Base.show(io::IO, mesh::UniformBZMesh)
    print(io, "Uniform BZ mesh (Cell = ", mesh.cell)
    print(io, ", origin = ", inv_lattice_vector(mesh) * mesh.origin)
    print(io, ", size = ", mesh.size)
    print(io, ", shift = ", mesh.shift)
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", mesh::UniformBZMesh)
    println(io, "Uniform BZ mesh:")
    if mesh.origin == 0
        println(io, "mesh type", "  Gamma-centered")
    else
        println(io, "mesh type", "  Monkhorst-Pack")
    end
    showfieldln(io, "origin = ", inv_lattice_vector(mesh) * mesh.origin)
    showfieldln(io, "size = ", mesh.size)
    showfieldln(io, "shift = ", mesh.shift)

    println(io)
    modelstr = sprint(show, "text/plain", mesh.cell)
    indent = " "^SHOWINDENTION
    print(io, indent, "Discretized " * replace(modelstr, "\n" => "\n" * indent))
end

################## Mesh Reduce ############################

# in spglib, grid_address runs from 1-ceil(N/2) to N-ceil(N/2)
# thus -1:2 for N=4 and -2:2 for N=5

function spglib_grid_address_to_index(mesh::AbstractUniformMesh{T,DIM}, ga) where {T,DIM}
    inds = ga[1:DIM] # if length(x)==3 but DIM==2, take first two
    fcoords = (inds .+ mesh.shift) ./ mesh.size #fractional coordinates as defined in spglib
    # shift fcoords, nomalize to [0, 1)
    fcoords = [(fcoords[i] < 0) ? (fcoords[i] + 1) : fcoords[i] for i in 1:DIM]
    x = lattice_vector(mesh) * fcoords
    return locate(mesh, x)
end

function MeshMaps.MeshMap(mesh::UniformBZMesh{T,DIM},
    # symmetry::Bool=true;
    is_time_reversal::Bool=true,
    tol_symmetry=PointSymmetry.SYMMETRY_TOLERANCE
) where {T,DIM}
    # Determine symmetry operations to use
    # if symmetry
    symmetries = Cells.default_symmetries(mesh.cell, tol_symmetry=tol_symmetry)
    # else
    #     symmetries = [one(PointSymmetry.SymOp)]
    # end
    @assert inv_lattice_vector(mesh) * mesh.origin ≈ -ones(DIM) / 2 "Uniform BZ mesh MeshMap only supports origin=[-1/2, -1/2,...], now got $(inv_lattice_vector(mesh) * mesh.origin)"

    kgrid_size = mesh.size
    kshift = mesh.shift
    @assert !isempty(symmetries)  # Identity has to be always present.

    _kgrid_size = ones(Int, 3)
    _kgrid_size[1:DIM] .= kgrid_size[1:DIM]
    is_shift = [0, 0, 0]
    is_shift[1:DIM] .= Int.(kshift * 2)
    # kcoords_mp = [mesh.mesh.inv_lattice * mesh[k] for k in 1:length(mesh)]
    # println("sym before: ", length(symmetries))
    # return kcoords_mp

    #TODO: It is probably makes sense to use the symmetry operations to reduce the number of kpoints
    # Filter those symmetry operations that preserve the MP grid
    # kcoords_mp = kgrid_monkhorst_pack(kgrid_size; kshift)
    # symmetries = symmetries_preserving_kgrid(symmetries, kcoords_mp)

    # Ws = [symop.W for symop in symmetries]
    # _, mapping, grid = PointSymmetry.spglib_get_stabilized_reciprocal_mesh(
    #     kgrid_size, Ws; is_shift, is_time_reversal=false
    # )

    # lat, atoms, pos, mag_moments = PointSymmetry.spglib_standardize_cell(mesh.cell, primitive=true)
    lat, atoms, pos, mag_moments = mesh.cell.lattice, mesh.cell.atoms, mesh.cell.positions, []

    lat, pos = PointSymmetry._make3D(lat, pos)

    cell, _ = PointSymmetry.spglib_cell(lat, atoms, pos, mag_moments)
    # println(cell)
    ngrid, mapping, _grid = PointSymmetry.get_ir_reciprocal_mesh(cell, _kgrid_size, is_shift;
        is_time_reversal=is_time_reversal, symprec=PointSymmetry.SYMMETRY_TOLERANCE)

    _grid = Int.(_grid)
    grid = Vector{Vector{Int}}()
    for i in 1:length(mesh)
        k = _grid[3*(i-1)+1:3*i]
        # println(k)
        push!(grid, k)
    end
    # println(grid)

    # if size = [4, 4, 4]
    # grid = 
    # [[ 0  0  0]
    #  [ 1  0  0]
    #  [ 2  0  0]
    #  [-1  0  0]
    #  [ 0  1  0]
    #  [ 1  1  0]
    #  [ 2  1  0]
    #  [-1  1  0]
    #  ....      ]
    @assert grid[2][1] - grid[1][1] == 1 "expect the first index to iterate first"

    # mapping .+= 1
    kidx_unique = unique(mapping)

    # Convert irreducible k-points to DFTK conventions
    kirreds = [(is_shift ./ 2 .+ grid[ik]) ./ _kgrid_size for ik in kidx_unique]
    # println(kirreds)

    # Find the indices of the corresponding reducible k-points in `grid`, which one of the
    # irreducible k-points in `kirreds` generates.
    k_all_reducible = [findall(isequal(elem), mapping) for elem in kidx_unique]
    # println(k_all_reducible)

    # Number of reducible k-points represented by the irreducible k-point `kirreds[ik]`
    n_equivalent_k = length.(k_all_reducible)
    @assert sum(n_equivalent_k) == prod(kgrid_size)
    kweights = n_equivalent_k / sum(n_equivalent_k)
    # println("kweights: ", kweights)

    # This loop checks for reducible k-points, which could not be mapped to any irreducible
    # k-point yet even though spglib claims this can be done.
    # This happens because spglib actually fails for some non-ideal lattices, resulting
    # in *wrong results* being returned. See the discussion in
    # https://github.com/spglib/spglib/issues/101
    for (iks_reducible, k) in zip(k_all_reducible, kirreds), ikred in iks_reducible
        kred = (is_shift ./ 2 .+ grid[ikred]) ./ _kgrid_size

        found_mapping = any(symmetries) do symop
            # If the difference between kred and W' * k == W^{-1} * k
            # is only integer in fractional reciprocal-space coordinates, then
            # kred and S' * k are equivalent k-points
            all(isinteger, kred - (symop.S * k))
        end

        if !found_mapping
            error("The reducible k-point $kred could not be generated from " *
                  "the irreducible kpoints. This points to a bug in spglib.")
        end
    end

    new_map = zeros(Int, length(mesh))
    for (i, m) in enumerate(mapping)
        grid_address = grid[i]
        mapped_grid_address = grid[m]
        kidx = spglib_grid_address_to_index(mesh, grid_address)
        mapped_kidx = spglib_grid_address_to_index(mesh, mapped_grid_address)
        new_map[kidx] = mapped_kidx
    end

    # return kirreds, mapping, k_all_reducible
    return MeshMaps.MeshMap(new_map)
end


# function _kcoords2ind(kcoord, kgrid_size, kshift)
#     # kidx = [Int(kvec * kgrid_size[d] - kshift[d]) + kgrid_size[d] - 1 for (d, kvec) in enumerate(kcoord)]
#     # inexact convert is not allowed with Int(), use floor(Int,) instead
#     kidx = [floor(Int, (kvec + 1 / 2) * kgrid_size[d] - kshift[d]) for (d, kvec) in enumerate(kcoord)]
#     kidx = [(kidx[d] + kgrid_size[d]) % kgrid_size[d] + 1 for d in 1:length(kidx)]
#     klinearidx = AbstractMeshes._inds2ind(tuple(kgrid_size...), tuple(kidx...))
#     return klinearidx
# end

# # WARNINING: Do not support Gamma_centered: origin=0  !!!
# # Monkhorst-Pack: origin=-1/2, consistent with VASP
# # to be consistent with DFTK: 
# #  - N is even, VASP is the same as DFTK: shift=0 will include Gamma point, shift=1/2 will not
# #  - N is odd, VASP is different as DFTK: shift=0 will not include Gamma point, shift=1/2 will
# function _reduced_uniform_meshmap(model::Cells.Cell{T,DIM}, symmetry::Bool=true;
#     kgrid_size::Vector{Int}, kshift::Bool=false,
#     tol_symmetry=PointSymmetry.SYMMETRY_TOLERANCE
# ) where {T,DIM}
#     # Determine symmetry operations to use
#     if symmetry
#         symmetries = default_symmetries(model, tol_symmetry=tol_symmetry)
#     else
#         symmetries = [one(PointSymmetry.SymOp)]
#     end
#     @assert !isempty(symmetries)  # Identity has to be always present.
#     _kgrid_size = ones(Int, 3)
#     _kgrid_size[1:DIM] = kgrid_size[1:DIM]
#     _kshift = kshift ? [1 // 2, 1 // 2, 1 // 2] : [0, 0, 0]

#     kcoords, kweights, symmetries = PointSymmetry.bzmesh_ir_wedge(_kgrid_size, symmetries; kshift=_kshift)
#     all_kcoords = PointSymmetry.unfold_kcoords(kcoords, symmetries)
#     Nk = reduce(*, kgrid_size)
#     @assert length(all_kcoords) == Nk
#     kindices = []
#     kmap = zeros(Int, Nk)
#     inv_kmap = Dict{Int,Vector{Int}}()
#     count = 0
#     for kpoint in kcoords
#         all_kpoint = PointSymmetry.unfold_kcoords([kpoint,], symmetries)
#         k0ind = _kcoords2ind(kpoint, kgrid_size, _kshift)
#         push!(kindices, k0ind)
#         inv_kmap[k0ind] = []
#         count += length(all_kpoint)
#         for k in all_kpoint
#             kind = _kcoords2ind(k, kgrid_size, _kshift)
#             kmap[kind] = k0ind
#             push!(inv_kmap[k0ind], kind)
#         end
#     end

#     @assert count == Nk

#     return MeshMaps.MeshMap(kindices, kmap, inv_kmap)
# end


include("PolarMeshes.jl")

end
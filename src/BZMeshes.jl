module BZMeshes

using ..StaticArrays
using ..LinearAlgebra

using ..AbstractMeshes
using ..AbstractMeshes: _inds2ind, _ind2inds
using ..Model
using ..BaseMesh
using ..BaseMesh: _indfloor
using ..BaseMesh.AbstractUniformMesh

export UniformBZMesh, DFTK_Monkhorst_Pack

"""
    struct UniformBZMesh{T, DIM} <: AbstractUniformMesh{T, DIM}

Uniformly distributed Brillouin zone mesh. Defined as a uniform mesh on 1st Brillouin zone
with Brillouin zone information stored in mesh.br::Brillouin. 

# Parameters:
- `T`: type of data
- `DIM`: dimension of the Brillouin zone

# Members:
- `br`: Brillouin zone information including lattice info, atom pos and allowed G vectors
- `origin`: origin of the uniform mesh, related to convention of 1st Brillouin zone. Commonly set to either (0,0,0) or such that (0,0,0) is at the center
- `size`: size of the uniform mesh. For Monkhorst-Pack mesh require even number.
- `shift`: k-shift of each mesh point. Set all to be zero for Gamma-centered and all to be 1//2 for M-P mesh
"""
# struct UniformBZMesh{T,DIM} <: AbstractMesh{T,DIM}
struct UniformBZMesh{T,DIM} <: AbstractUniformMesh{T,DIM}
    br::Brillouin{T,DIM}

    origin::SVector{DIM,T}
    size::NTuple{DIM,Int}
    shift::SVector{DIM,Rational}
end

# default shift is 1/2, result in Monkhorst-Pack mesh
# with shift = 0, result in Gamma-centered
# can also customize with shift::SVector by calling default constructor
# \Gamma=(0,0,0) is at center by default, can be set at corner by setting origin to it
"""
    function UniformBZMesh(; br::Brillouin, origin, size, shift)

customized constructor for UniformBZMesh. The parameters origin and shift is provided to customize
the mesh as Gamma-centered or M-P mesh. 

# Parameters:
- `br`: Brillouin zone info
- `origin`: a number indicating shift of origin. 
    the actuall origin becomes origin*(b1+b2+b3)
    default value origin=-1/2 takes (0,0,0) to center of 1st BZ, origin=0 makes mesh[1,1,1]=(0,0,0)+shift
- `size`: size of the mesh
- `shift`: additional k-shift for mesh points. 
    actuall shift is shift*(b1/N1+b2/N2+b3/N3)
    for even N, shift=1/2 avoids high symmetry points while preserve symmetry.
"""
UniformBZMesh(;
    br::Brillouin{T,DIM},
    origin::Real=-1 // 2,
    size,
    shift::Real=1 // 2) where {T,DIM} = UniformBZMesh{T,DIM}(
    br, origin .* ones(T, DIM), size, shift .* ones(Int, DIM)
)

function DFTK_Monkhorst_Pack(;
    br::Brillouin{T,DIM},
    size,
    shift::AbstractVector) where {T,DIM}
    kshift = [(iseven(size[i]) ? shift[i] : shift[i] + 1 // 2) for i in 1:DIM]
    return UniformBZMesh{T,DIM}(
        br,
        -1 // 2 .* ones(T, DIM), # origin
        size,
        kshift
    )
end

function Monkhorst_Pack(;
    # Gamma_centered: origin=0, 
    # Monkhorst-Pack: origin=-1/2, consistent with VASP
    # to be consistent with DFTK: 
    #  - N is even, VASP is the same as DFTK: shift=0 will include Gamma point, shift=1/2 will not
    #  - N is odd, VASP is different as DFTK: shift=0 will not include Gamma point, shift=1/2 will
    br::Brillouin{T,DIM},
    size,
    shift::AbstractVector=[0, 0, 0]
) where {T,DIM}
    # kshift = [(iseven(size[i]) ? shift[i] : shift[i] + 1 // 2) for i in 1:DIM]
    return UniformBZMesh{T,DIM}(
        br,
        -ones(T, DIM) / 2,
        tuple(size...),
        shift=shift
    )
end

lattice_vector(mesh::UniformBZMesh) = mesh.br.recip_lattice
inv_lattice_vector(mesh::UniformBZMesh) = mesh.inv_recip_lattice
cell_volume(mesh::UniformBZMesh) = mesh.br.recip_cell_volume

function Base.show(io::IO, mesh::UniformBZMesh)
    println("UniformBZMesh with $(length(mesh)) mesh points")
end

end
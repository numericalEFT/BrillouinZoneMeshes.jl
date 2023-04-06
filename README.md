# BrillouinZoneMeshes

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://numericalEFT.github.io/BrillouinZoneMeshes.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://numericalEFT.github.io/BrillouinZoneMeshes.jl/dev)
[![Build Status](https://github.com/numericalEFT/BrillouinZoneMeshes.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/numericalEFT/BrillouinZoneMeshes.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/numericalEFT/BrillouinZoneMeshes.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/numericalEFT/BrillouinZoneMeshes.jl)

# BrillouinZoneMeshes

Documentation for [BrillouinZoneMeshes](https://github.com/numericalEFT/BrillouinZoneMeshes.jl).

This package provides general-purpose multi-dimensional meshes for numerical 
representation of continuous functions and specialized meshes 
for functions on Brillouin Zones. 

## Getting started

Setup with:

```julia
    DIM = 2;
    # square lattice
    lattice = Matrix([1.0 0; 0 1]');
    # create Brillouin zone
    br = BZMeshes.Cell(lattice=lattice);
    # uniform mesh
    umesh = BZMeshes.UniformBZMesh(cell=br, size=(4, 4));
    # symmetry reduce map
    mm = MeshMaps.MeshMap(umesh);
    # reduced mesh
    rmesh = ReducedBZMesh(umesh, mm);
```

and here are some examples of usage:

```julia
julia> using BrillouinZoneMeshes; DIM = 2; lattice = Matrix([1.0 0; 0 1]'); br = BZMeshes.Cell(lattice=lattice); umesh = BZMeshes.UniformBZMesh(cell=br, size=(4, 4)); mm = MeshMaps.MeshMap(umesh); rmesh = ReducedBZMesh(umesh, mm);

julia> length(umesh), length(rmesh)
(16, 3)

julia> AbstractMeshes.locate(rmesh, [1,1])
3

julia> rmesh[3]
2-element StaticArraysCore.SVector{2, Float64} with indices SOneTo(2):
 0.7853981633974483
 0.7853981633974483

julia> AbstractMeshes.volume(rmesh) / 4π^2
1.0

julia> data = ones(3)
3-element Vector{Float64}:
 1.0
 1.0
 1.0

julia> AbstractMeshes.integrate(data, rmesh)
39.47841760435743

julia> AbstractMeshes.interp(data, rmesh, [0.3,-0.2])
1.0
```

## General

Various mesh grids for different purposes are defined as concrete types 
derived from `AbstractMeshes.AbstractMesh`. All of them are supposed to 
behave as `AbstractArray` with elements being `SVector` representing the
mesh points in Cartesian coordinates. 

In addition to the interface of `AbstractMeshes.AbstractMesh`, four useful
methods are defined: `locate`, `volume`, `interp`, and `integrate`. 
* `locate(mesh, x)` finds the mesh point nearest to x
* `volume(mesh, i)` gives the volume represented by mesh point `mesh[i]`
* `interp(data, mesh, x)` gives the interpolation of `data` on `mesh` at    point `x`
* `integrate(data, mesh)` compute integration of `data` on `mesh`

If it is known that some of the mesh points are guaranteed to have the same 
data value, it's possible to define a `MeshMap` to reveal this fact and 
create a `ReducedBZMesh` to save storage space.

## Brillouin zone

The information of Brillouin zone is stored in `Cells.Cell`. 
Including lattice vector, reciprocal lattice vector and their inverse;
volume of unit cell and reciprocal unit cell; G vectors for extended 
Brillouin zone and symmetries.

## Uniform Meshes

Uniform meshes are defined as uniformly distributed meshes on a 
parallellogram area described by an origin and a set of lattice vectors.
The simplest one is `BaseMesh.UMesh`, while `BZMeshes.UniformBZMesh` 
containes additional information about the Brillouin zone stored in
its `cell` field. 

Uniform meshes are conventionally used in various _ab initio_ calculations. 
In this package various frequently used meshes, such as Gamma-centered and 
Monkhorst-Pack meshes, could be generated via `BZMeshes.UniformBZMesh` with 
different parameters. The default parameter of the constructor of 
`BZMeshes.UniformBZMesh` generates Gamma-centered mesh, while two constructors
for M-P mesh, `Monkhorst_Pack` and `DFTK_Monkhorst_Pack`, follow conventions
from VASP and DFTK respectively.

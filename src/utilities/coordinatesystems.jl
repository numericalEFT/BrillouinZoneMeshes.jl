module Coordinates

using ..StaticArrays
using ..LinearAlgebra
# copied from CoordinateTransformations.jl, modify some conventions for our usage

# Core methods
export compose, ∘, transform_deriv, transform_deriv_params, recenter
export Transformation, IdentityTransformation

# 2D coordinate systems and their transformations
export Polar
export PolarFromCartesian, CartesianFromPolar

# 3D coordinate systems and their transformations
export Spherical, Cylindrical
export SphericalFromCartesian, CartesianFromSpherical,
    CylindricalFromCartesian, CartesianFromCylindrical,
    CylindricalFromSpherical, SphericalFromCylindrical

#####################################
### Interface for transformations ###
#####################################

# The external interface consists of transform(), Base.inv(),  compose() (or ∘),
# transform_deriv() and transform_deriv_params()

"""
The `Transformation` supertype defines a simple interface for performing
transformations. Subtypes should be able to apply a coordinate system
transformation on the correct data types by overloading the call method, and
usually would have the corresponding inverse transformation defined by `Base.inv()`.
Efficient compositions can optionally be defined by `compose()` (equivalently `∘`).
"""
abstract type Transformation end

"""
The `IdentityTransformation` is a singleton `Transformation` that returns the
input unchanged, similar to `identity`.
"""
struct IdentityTransformation <: Transformation end

@inline (::IdentityTransformation)(x) = x

"""
A `ComposedTransformation` simply executes two transformations successively, and
is the fallback output type of `compose()`.
"""
struct ComposedTransformation{T1<:Transformation,T2<:Transformation} <: Transformation
    t1::T1
    t2::T2
end

Base.show(io::IO, trans::ComposedTransformation) = print(io, "($(trans.t1) ∘ $(trans.t2))")

@inline function (trans::ComposedTransformation)(x)
    trans.t1(trans.t2(x))
end

function Base.:(==)(trans1::ComposedTransformation, trans2::ComposedTransformation)
    (trans1.t1 == trans2.t1) && (trans1.t2 == trans2.t2)
end

function Base.isapprox(trans1::ComposedTransformation, trans2::ComposedTransformation; kwargs...)
    isapprox(trans1.t1, trans2.t1; kwargs...) && isapprox(trans1.t2, trans2.t2; kwargs...)
end

const compose = ∘

"""
    compose(trans1, trans2)
    trans1 ∘ trans2
Take two transformations and create a new transformation that is equivalent to
successively applying `trans2` to the coordinate, and then `trans1`. By default
will create a `ComposedTransformation`, however this method can be overloaded
for efficiency (e.g. two affine transformations naturally compose to a single
affine transformation).
"""
function compose(trans1::Transformation, trans2::Transformation)
    ComposedTransformation(trans1, trans2)
end

compose(trans::IdentityTransformation, ::IdentityTransformation) = trans
compose(::IdentityTransformation, trans::Transformation) = trans
compose(trans::Transformation, ::IdentityTransformation) = trans


"""
    inv(trans::Transformation)
Returns the inverse (or reverse) of the transformation `trans`
"""
function Base.inv(trans::Transformation)
    error("Inverse transformation for $(typeof(trans)) has not been defined.")
end

Base.inv(trans::ComposedTransformation) = inv(trans.t2) ∘ inv(trans.t1)
Base.inv(trans::IdentityTransformation) = trans

"""
    recenter(trans::Union{AbstractMatrix,Transformation}, origin::Union{AbstractVector, Tuple}) -> ctrans
Return a new transformation `ctrans` such that point `origin` serves
as the origin-of-coordinates for `trans`. Translation by `±origin`
occurs both before and after applying `trans`, so that if `trans` is
linear we have
    ctrans(origin) == origin
As a consequence, `recenter` only makes sense if the output space of
`trans` is isomorϕc with the input space.
For example, if `trans` is a rotation matrix, then `ctrans` rotates
space around `origin`.
"""
function recenter(trans::Transformation, origin::AbstractVector)
    Translation(origin) ∘ trans ∘ Translation(-origin)
end
recenter(trans::Transformation, origin::Tuple) = recenter(trans, SVector(origin))


"""
    transform_deriv(trans::Transformation, x)
A matrix describing how differentials on the parameters of `x` flow through to
the output of transformation `trans`.
"""
transform_deriv(trans::Transformation, x) = error("Differential matrix of transform $trans with input $x not defined")

transform_deriv(::IdentityTransformation, x) = I

function transform_deriv(trans::ComposedTransformation, x)
    x2 = trans.t2(x)
    m1 = transform_deriv(trans.t1, x2)
    m2 = transform_deriv(trans.t2, x)
    return m1 * m2
end


"""
    transform_deriv_params(trans::Transformation, x)
A matrix describing how differentials on the parameters of `trans` flow through
to the output of transformation `trans` given input `x`.
"""
transform_deriv_params(trans::Transformation, x) = error("Differential matrix of parameters of transform $trans with input $x not defined")

transform_deriv_params(::IdentityTransformation, x) = error("IdentityTransformation has no parameters")

function transform_deriv_params(trans::ComposedTransformation, x)
    x2 = trans.t2(x)
    m1 = transform_deriv(trans.t1, x2)
    p2 = transform_deriv_params(trans.t2, x)
    p1 = transform_deriv_params(trans.t1, x2)
    return hcat(p1, m1 * p2)
end

#############################
### 2D Coordinate systems ###
#############################
"""
`Polar{T,A}(r::T, ϕ::A)` - 2D polar coordinates
"""
struct Polar{T,A}
    r::T
    ϕ::A

    Polar{T,A}(r, ϕ) where {T,A} = new(r, ϕ)
end

function Polar(r, ϕ)
    r2, ϕ2 = promote(r, ϕ)

    return Polar{typeof(r2),typeof(ϕ2)}(r2, ϕ2)
end

Base.show(io::IO, x::Polar) = print(io, "Polar(r=$(x.r), ϕ=$(x.ϕ) rad)")
Base.isapprox(p1::Polar, p2::Polar; kwargs...) = isapprox(p1.r, p2.r; kwargs...) && isapprox(p1.ϕ, p2.ϕ; kwargs...)

"`PolarFromCartesian()` - transformation from `AbstractVector` of length 2 to `Polar` type"
struct PolarFromCartesian <: Transformation end
"`CartesianFromPolar()` - transformation from `Polar` type to `SVector{2}` type"
struct CartesianFromPolar <: Transformation end

Base.show(io::IO, trans::PolarFromCartesian) = print(io, "PolarFromCartesian()")
Base.show(io::IO, trans::CartesianFromPolar) = print(io, "CartesianFromPolar()")

function (::PolarFromCartesian)(x::AbstractVector)
    length(x) == 2 || error("Polar transform takes a 2D coordinate")

    Polar(hypot(x[1], x[2]), atan(x[2], x[1]))
end

function transform_deriv(::PolarFromCartesian, x::AbstractVector)
    length(x) == 2 || error("Polar transform takes a 2D coordinate")

    r = hypot(x[1], x[2])
    f = x[2] / x[1]
    c = one(eltype(x)) / (x[1] * (one(eltype(x)) + f * f))
    @SMatrix [x[1]/r x[2]/r
        -f*c c]
end
transform_deriv_params(::PolarFromCartesian, x::AbstractVector) = error("PolarFromCartesian has no parameters")

function (::CartesianFromPolar)(x::Polar)
    s, c = sincos(x.ϕ)
    SVector(x.r * c, x.r * s)
end
function transform_deriv(::CartesianFromPolar, x::Polar)
    sϕ, cϕ = sincos(x.ϕ)
    @SMatrix [cϕ -x.r*sϕ
        sϕ x.r*cϕ]
end
transform_deriv_params(::CartesianFromPolar, x::Polar) = error("CartesianFromPolar has no parameters")

Base.inv(::PolarFromCartesian) = CartesianFromPolar()
Base.inv(::CartesianFromPolar) = PolarFromCartesian()

compose(::PolarFromCartesian, ::CartesianFromPolar) = IdentityTransformation()
compose(::CartesianFromPolar, ::PolarFromCartesian) = IdentityTransformation()

# For convenience
Base.convert(::Type{Polar}, v::AbstractVector) = PolarFromCartesian()(v)
@inline Base.convert(::Type{V}, p::Polar) where {V<:AbstractVector} = convert(V, CartesianFromPolar()(p))
@inline Base.convert(::Type{V}, p::Polar) where {V<:StaticVector} = convert(V, CartesianFromPolar()(p))


#############################
### 3D Coordinate Systems ###
#############################
"""
Spherical(r, θ, ϕ) - 3D spherical coordinates

There are many Spherical coordinate conventions and this library uses a somewhat exotic one.
Given a vector `v` with Cartesian coordinates `xyz`, let `v_xy = [x,y,0]` be the orthogonal projection of `v` on the `xy` plane.

* `r` is the radius. It is given by `norm(v, 2)`.
* `θ` is the latitude. It is the angle from `v_xy` to `v`.
* `ϕ` is the azimuth. It is the angle from the x-axis to `v_xy`

```jldoctest
julia> using CoordinateTransformations

julia> v = randn(3);

julia> sph = SphericalFromCartesian()(v);

julia> r = sph.r; ϕ=sph.ϕ; θ=sph.θ;

julia> v ≈ [r * cos(ϕ) * cos(θ), r * sin(ϕ) * cos(θ), r*sin(θ)]
true
"""
struct Spherical{T,A}
    r::T
    θ::A
    ϕ::A

    Spherical{T,A}(r, θ, ϕ) where {T,A} = new(r, θ, ϕ)
end

function Spherical(r, θ, ϕ)
    r2, θ2, ϕ2 = promote(r, θ, ϕ)

    return Spherical{typeof(r2),typeof(ϕ2)}(r2, θ2, ϕ2)
end

Base.show(io::IO, x::Spherical) = print(io, "Spherical(r=$(x.r), θ=$(x.θ) rad, ϕ=$(x.ϕ) rad)")
Base.isapprox(p1::Spherical, p2::Spherical; kwargs...) = isapprox(p1.r, p2.r; kwargs...) && isapprox(p1.ϕ, p2.ϕ; kwargs...) && isapprox(p1.θ, p2.θ; kwargs...)

"""
Cylindrical(r, ϕ, z) - 3D cylindrical coordinates
"""
struct Cylindrical{T,A}
    r::T
    ϕ::A
    z::T

    Cylindrical{T,A}(r, ϕ, z) where {T,A} = new(r, ϕ, z)
end

function Cylindrical(r, ϕ, z)
    r2, ϕ2, z2 = promote(r, ϕ, z)

    return Cylindrical{typeof(r2),typeof(ϕ2)}(r2, ϕ2, z2)
end

Base.show(io::IO, x::Cylindrical) = print(io, "Cylindrical(r=$(x.r), ϕ=$(x.ϕ) rad, z=$(x.z))")
Base.isapprox(p1::Cylindrical, p2::Cylindrical; kwargs...) = isapprox(p1.r, p2.r; kwargs...) && isapprox(p1.ϕ, p2.ϕ; kwargs...) && isapprox(p1.z, p2.z; kwargs...)

"`SphericalFromCartesian()` - transformation from 3D point to `Spherical` type"
struct SphericalFromCartesian <: Transformation end
"`CartesianFromSpherical()` - transformation from `Spherical` type to `SVector{3}` type"
struct CartesianFromSpherical <: Transformation end
"`CylindricalFromCartesian()` - transformation from 3D point to `Cylindrical` type"
struct CylindricalFromCartesian <: Transformation end
"`CartesianFromCylindrical()` - transformation from `Cylindrical` type to `SVector{3}` type"
struct CartesianFromCylindrical <: Transformation end
"`CylindricalFromSpherical()` - transformation from `Spherical` type to `Cylindrical` type"
struct CylindricalFromSpherical <: Transformation end
"`SphericalFromCylindrical()` - transformation from `Cylindrical` type to `Spherical` type"
struct SphericalFromCylindrical <: Transformation end

Base.show(io::IO, trans::SphericalFromCartesian) = print(io, "SphericalFromCartesian()")
Base.show(io::IO, trans::CartesianFromSpherical) = print(io, "CartesianFromSpherical()")
Base.show(io::IO, trans::CylindricalFromCartesian) = print(io, "CylindricalFromCartesian()")
Base.show(io::IO, trans::CartesianFromCylindrical) = print(io, "CartesianFromCylindrical()")
Base.show(io::IO, trans::CylindricalFromSpherical) = print(io, "CylindricalFromSpherical()")
Base.show(io::IO, trans::SphericalFromCylindrical) = print(io, "SphericalFromCylindrical()")

# Cartesian <-> Spherical
function (::SphericalFromCartesian)(x::AbstractVector)
    length(x) == 3 || error("Spherical transform takes a 3D coordinate")

    Spherical(hypot(x[1], x[2], x[3]), atan(x[3], hypot(x[1], x[2])), atan(x[2], x[1]))
end

function (::CartesianFromSpherical)(x::Spherical)
    sϕ, cϕ = sincos(x.ϕ)
    sθ, cθ = sincos(x.θ)
    SVector(x.r * cϕ * cθ, x.r * sϕ * cθ, x.r * sθ)
end

# Cartesian <-> Cylindrical
function (::CylindricalFromCartesian)(x::AbstractVector)
    length(x) == 3 || error("Cylindrical transform takes a 3D coordinate")

    Cylindrical(hypot(x[1], x[2]), atan(x[2], x[1]), x[3])
end

function transform_deriv(::CylindricalFromCartesian, x::AbstractVector)
    length(x) == 3 || error("Cylindrical transform takes a 3D coordinate")
    T = eltype(x)

    r = hypot(x[1], x[2])
    f = x[2] / x[1]
    c = one(T) / (x[1] * (one(T) + f * f))
    @SMatrix [x[1]/r x[2]/r zero(T)
        -f*c c zero(T)
        zero(T) zero(T) one(T)]
end
transform_deriv_params(::CylindricalFromCartesian, x::AbstractVector) = error("CylindricalFromCartesian has no parameters")

function (::CartesianFromCylindrical)(x::Cylindrical)
    sϕ, cϕ = sincos(x.ϕ)
    SVector(x.r * cϕ, x.r * sϕ, x.z)
end
function transform_deriv(::CartesianFromCylindrical, x::Cylindrical{T}) where {T}
    sϕ, cϕ = sincos(x.ϕ)
    @SMatrix [cϕ -x.r*sϕ zero(T)
        sϕ x.r*cϕ zero(T)
        zero(T) zero(T) one(T)]
end
transform_deriv_params(::CartesianFromPolar, x::Cylindrical) = error("CartesianFromCylindrical has no parameters")

function (::CylindricalFromSpherical)(x::Spherical)
    sθ, cθ = sincos(x.θ)
    Cylindrical(x.r * cθ, x.ϕ, x.r * sθ)
end
function transform_deriv(::CylindricalFromSpherical, x::Spherical)
    M1 = transform_deriv(CylindricalFromCartesian(), CartesianFromSpherical()(x))
    M2 = transform_deriv(CartesianFromSpherical(), x)
    return M1 * M2
end
transform_deriv_params(::CylindricalFromSpherical, x::Spherical) = error("CylindricalFromSpherical has no parameters")

function (::SphericalFromCylindrical)(x::Cylindrical)
    Spherical(hypot(x.r, x.z), x.ϕ, atan(x.z, x.r))
end
function transform_deriv(::SphericalFromCylindrical, x::Cylindrical)
    M1 = transform_deriv(SphericalFromCartesian(), CartesianFromCylindrical()(x))
    M2 = transform_deriv(CartesianFromCylindrical(), x)
    return M1 * M2
end
transform_deriv_params(::SphericalFromCylindrical, x::Cylindrical) = error("SphericalFromCylindrical has no parameters")

Base.inv(::SphericalFromCartesian) = CartesianFromSpherical()
Base.inv(::CartesianFromSpherical) = SphericalFromCartesian()
Base.inv(::CylindricalFromCartesian) = CartesianFromCylindrical()
Base.inv(::CartesianFromCylindrical) = CylindricalFromCartesian()
Base.inv(::CylindricalFromSpherical) = SphericalFromCylindrical()
Base.inv(::SphericalFromCylindrical) = CylindricalFromSpherical()

# Inverse composition
compose(::SphericalFromCartesian, ::CartesianFromSpherical) = IdentityTransformation()
compose(::CartesianFromSpherical, ::SphericalFromCartesian) = IdentityTransformation()
compose(::CylindricalFromCartesian, ::CartesianFromCylindrical) = IdentityTransformation()
compose(::CartesianFromCylindrical, ::CylindricalFromCartesian) = IdentityTransformation()
compose(::CylindricalFromSpherical, ::SphericalFromCylindrical) = IdentityTransformation()
compose(::SphericalFromCylindrical, ::CylindricalFromSpherical) = IdentityTransformation()

# Cyclic compositions
compose(::SphericalFromCartesian, ::CartesianFromCylindrical) = SphericalFromCylindrical()
compose(::CartesianFromSpherical, ::SphericalFromCylindrical) = CartesianFromCylindrical()
compose(::CylindricalFromCartesian, ::CartesianFromSpherical) = CylindricalFromSpherical()
compose(::CartesianFromCylindrical, ::CylindricalFromSpherical) = CartesianFromSpherical()
compose(::CylindricalFromSpherical, ::SphericalFromCartesian) = CylindricalFromCartesian()
compose(::SphericalFromCylindrical, ::CylindricalFromCartesian) = SphericalFromCartesian()

# For convenience
Base.convert(::Type{Spherical}, v::AbstractVector) = SphericalFromCartesian()(v)
Base.convert(::Type{Cylindrical}, v::AbstractVector) = CylindricalFromCartesian()(v)

Base.convert(::Type{V}, s::Spherical) where {V<:AbstractVector} = convert(V, CartesianFromSpherical()(s))
Base.convert(::Type{V}, c::Cylindrical) where {V<:AbstractVector} = convert(V, CartesianFromCylindrical()(c))
Base.convert(::Type{V}, s::Spherical) where {V<:StaticVector} = convert(V, CartesianFromSpherical()(s))
Base.convert(::Type{V}, c::Cylindrical) where {V<:StaticVector} = convert(V, CartesianFromCylindrical()(c))

Base.convert(::Type{Spherical}, c::Cylindrical) = SphericalFromCylindrical()(c)
Base.convert(::Type{Cylindrical}, s::Spherical) = CylindricalFromSpherical()(s)

end
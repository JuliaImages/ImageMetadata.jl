module ImageMetadata

# The order here is designed to avoid an ambiguity warning in convert,
# see the top of ImageAxes
using ImageAxes
using ImageCore, Colors, FixedPointNumbers
using ColorVectorSpace   # for overriding math operations with Gray/RGB
import AxisArrays

import Base: +, -, *, /
import Base: permutedims

const ViewIndex = Union{Base.ViewIndex, Colon}

export
    # types
    ImageMeta,

    # functions
    copyproperties,
    data,
    properties,
    shareproperties,
    spatialproperties

#### types and constructors ####

# Concrete types
"""
`ImageMeta` is an AbstractArray that can have metadata, stored in a dictionary.

Construct an image with `ImageMeta(A, props)` (for a properties dictionary
`props`), or with `ImageMeta(A, prop1=val1, prop2=val2, ...)`.
"""
mutable struct ImageMeta{T,N,A<:AbstractArray,P<:AbstractDict{Symbol,Any}} <: AbstractArray{T,N}
    data::A
    properties::P

    function ImageMeta{T,N,A,P}(data::AbstractArray, properties::P) where {T,N,A,P}
        new{T,N,A,P}(data, properties)
    end
end

function ImageMeta(data::AbstractArray{T,N}, props::AbstractDict{Symbol,Any}) where {T,N}
    return ImageMeta{T,N,typeof(data),typeof(props)}(data,props)
end
function ImageMeta(data::AbstractArray{T,N}, props::Dict{Symbol,Any}) where {T,N}
    return ImageMeta{T,N,typeof(data),typeof(props)}(data,props)
end
function ImageMeta(data::AbstractArray, props::Dict{Symbol})
    return ImageMeta(data, convert(Dict{Symbol,Any}, props))
end
function ImageMeta(data::AbstractArray; kwargs...)
    return ImageMeta(data, kwargs2dict(kwargs))
end
ImageMeta(data::AbstractArray, props::AbstractDict) = throw(ArgumentError("properties must be an AbstractDict{Symbol,Any}"))

const ImageMetaArray{T,N,A<:Array} = ImageMeta{T,N,A}
const ImageMetaAxis{T,N,A<:AxisArray} = ImageMeta{T,N,A}

Base.size(A::ImageMeta) = size(data(A))
Base.size(A::ImageMetaAxis, Ax::Axis) = size(data(A), Ax)
Base.size(A::ImageMetaAxis, ::Type{Ax}) where {Ax<:Axis} = size(data(A), Ax)
Base.axes(A::ImageMeta) = axes(data(A))
Base.axes(A::ImageMetaAxis, Ax::Axis) = axes(data(A), Ax)
Base.axes(A::ImageMetaAxis, ::Type{Ax}) where {Ax<:Axis} = axes(data(A), Ax)

datatype(::Type{ImageMeta{T,N,A,P}}) where {T,N,A<:AbstractArray,P} = A

Base.IndexStyle(::Type{M}) where {M<:ImageMeta} = IndexStyle(datatype(M))

AxisArrays.HasAxes(A::ImageMetaAxis) = AxisArrays.HasAxes{true}()

# getindex and setindex!
for AType in (ImageMeta, ImageMetaAxis)
    @eval begin
        @inline function Base.getindex(img::$AType{T,1}, i::Int) where T
            @boundscheck checkbounds(data(img), i)
            @inbounds ret = data(img)[i]
            ret
        end
        @inline function Base.getindex(img::$AType, i::Int)
            @boundscheck checkbounds(data(img), i)
            @inbounds ret = data(img)[i]
            ret
        end
        @inline function Base.getindex(img::$AType{T,N}, I::Vararg{Int,N}) where {T,N}
            @boundscheck checkbounds(data(img), I...)
            @inbounds ret = data(img)[I...]
            ret
        end

        @inline function Base.setindex!(img::$AType{T,1}, val, i::Int) where T
            @boundscheck checkbounds(data(img), i)
            @inbounds data(img)[i] = val
            val
        end
        @inline function Base.setindex!(img::$AType, val, i::Int)
            @boundscheck checkbounds(data(img), i)
            @inbounds data(img)[i] = val
            val
        end
        @inline function Base.setindex!(img::$AType{T,N}, val, I::Vararg{Int,N}) where {T,N}
            @boundscheck checkbounds(data(img), I...)
            @inbounds data(img)[I...] = val
            val
        end
    end
end

@inline function Base.getindex(img::ImageMetaAxis, ax::Axis, I...)
    result = data(img)[ax, I...]
    maybe_wrap(img, result)
end
@inline function Base.getindex(img::ImageMetaAxis, i::Union{Integer,AbstractVector,Colon}, I...)
    result = data(img)[i, I...]
    maybe_wrap(img, result)
end
maybe_wrap(img::ImageMeta{T}, result::T) where T = result
maybe_wrap(img::ImageMeta{T}, result::AbstractArray{T}) where T = copyproperties(img, result)


@inline function Base.setindex!(img::ImageMetaAxis, val, ax::Axis, I...)
    setindex!(data(img), val, ax, I...)
end
@inline function Base.setindex!(img::ImageMetaAxis, val, i::Union{Integer,AbstractVector,Colon}, I...)
    setindex!(data(img), val, i, I...)
end

Base.view(img::ImageMeta, ax::Axis, I...) = shareproperties(img, view(data(img), ax, I...))
Base.view(img::ImageMeta{T,N}, I::Vararg{ViewIndex,N}) where {T,N} = shareproperties(img, view(data(img), I...))
Base.view(img::ImageMeta, i::ViewIndex) = shareproperties(img, view(data(img), i))
Base.view(img::ImageMeta, I::Vararg{ViewIndex,N}) where {N} = shareproperties(img, view(data(img), I...))

function Base.getproperty(img::ImageMeta, propname::Symbol)
    # TODO remove these once deprecations are done
    if propname === :data
        return getproperty_data(img)
    elseif propname === :properties
        return getproperty_properties(img)
    else
        return properties(img)[propname]
    end
end

function Base.setproperty!(img::ImageMeta, propname::Symbol, X)
    return setindex!(properties(img), X, propname )
end

Base.propertynames(img::ImageMeta) = (keys(properties(img))...,)

Base.copy(img::ImageMeta) = ImageMeta(copy(data(img)), deepcopy(properties(img)))

Base.convert(::Type{ImageMeta}, A::ImageMeta) = A
Base.convert(::Type{ImageMeta}, A::AbstractArray) = ImageMeta(A)
Base.convert(::Type{ImageMeta{T}}, A::ImageMeta{T}) where {T} = A
Base.convert(::Type{ImageMeta{T}}, A::ImageMeta) where {T} = shareproperties(A, convert(Array{T}, data(A)))
Base.convert(::Type{ImageMeta{T}}, A::AbstractArray) where {T} = ImageMeta(convert(Array{T}, A))

# copy properties
function Base.copy!(imgdest::ImageMeta, imgsrc::ImageMeta, prop1::Symbol, props::Symbol...)
    setproperty!(imgdest, prop1, deepcopy(getproperty(imgsrc, prop1)))
    for p in props
        setproperty!(imgdest, p, deepcopy(getproperty(imgsrc, p)))
    end
    return imgdest
end

# similar
Base.similar(img::ImageMeta, ::Type{T}, shape::Dims) where {T} = ImageMeta(similar(data(img), T, shape), deepcopy(properties(img)))
Base.similar(img::ImageMetaAxis, ::Type{T}) where {T} = ImageMeta(similar(data(img), T), deepcopy(properties(img)))

"""
    copyproperties(img::ImageMeta, data) -> imgnew

Create a new "image," copying the properties dictionary of `img` but
using the data of the AbstractArray `data`. Note that changing the
properties of `imgnew` does not affect the properties of `img`.

See also: [`shareproperties`](@ref).
"""
copyproperties(img::ImageMeta, data::AbstractArray) =
    ImageMeta(data, deepcopy(properties(img)))

"""
    shareproperties(img::ImageMeta, data) -> imgnew

Create a new "image," reusing the properties dictionary of `img` but
using the data of the AbstractArray `data`. The two images have
synchronized properties; modifying one also affects the other.

See also: [`copyproperties`](@ref).
"""
shareproperties(img::ImageMeta, data::AbstractArray) = ImageMeta(data, properties(img))

# Delete a property!
Base.delete!(img::ImageMeta, propname::Symbol) = delete!(properties(img), propname)


# Iteration
# Defer to the array object in case it has special iteration defined
Base.iterate(img::ImageMeta) = Base.iterate(data(img))
Base.iterate(img::ImageMeta, s) = Base.iterate(data(img), s)

# Show
const emptyset = Set()
function showim(io::IO, img::ImageMeta)
    IT = typeof(img)
    print(io, eltype(img).name.name, " ImageMeta with:\n  data: ", summary(data(img)), "\n  properties:")
    showdictlines(io, properties(img), get(img, :suppress, emptyset))
end
Base.show(io::IO, img::ImageMeta) = showim(io, img)
Base.show(io::IO, ::MIME"text/plain", img::ImageMeta) = showim(io, img)

function Base.reinterpret(::Type{T}, img::ImageMeta) where {T}
    shareproperties(img, reinterpret(T, data(img)))
end


"""
    data(img::ImageMeta) -> array

Extract the data from `img`, omitting the properties
dictionary. `array` shares storage with `img`, so changes to one
affect the other.

See also: [`properties`](@ref).
"""
ImageAxes.data(img::ImageMeta) = getfield(img, :data)  # fixme when deprecation is removed from ImageCore

function ImageCore.permuteddimsview(A::ImageMeta, perm)
    ip = sortperm([perm...][[coords_spatial(A)...]])  # the inverse spatial permutation
    permutedims_props!(copyproperties(A, permuteddimsview(data(A), perm)), ip)
end
ImageCore.channelview(A::ImageMeta) = shareproperties(A, channelview(data(A)))
ImageCore.colorview(::Type{C}, A::ImageMeta{T,N}) where {C<:Colorant,T,N} = shareproperties(A, colorview(C, data(A)))
ImageCore.colorview(::Type{ARGB32}, A::ImageMeta{T,N}) where {T,N} = shareproperties(A, colorview(ARGB32, data(A)))
ImageCore.rawview(A::ImageMeta{T}) where {T<:Real} = shareproperties(A, rawview(data(A)))
ImageCore.normedview(::Type{T}, A::ImageMeta{S}) where {T<:FixedPoint,S<:Unsigned} = shareproperties(A, normedview(T, data(A)))

# AxisArrays functions
AxisArrays.axes(img::ImageMetaAxis) = AxisArrays.axes(data(img))
AxisArrays.axes(img::ImageMetaAxis, d::Int) = AxisArrays.axes(data(img), d)
AxisArrays.axes(img::ImageMetaAxis, Ax::Axis) = AxisArrays.axes(data(img), Ax)
AxisArrays.axes(img::ImageMetaAxis, ::Type{Ax}) where {Ax<:Axis} = AxisArrays.axes(data(img), Ax)
AxisArrays.axisdim(img::ImageMetaAxis, ax) = axisdim(data(img), ax)
AxisArrays.axisnames(img::ImageMetaAxis) = axisnames(data(img))
AxisArrays.axisvalues(img::ImageMetaAxis) = axisvalues(data(img))

#### Properties ####

"""
    properties(imgmeta) -> props

Extract the properties dictionary `props` for `imgmeta`. `props`
shares storage with `img`, so changes to one affect the other.

See also: [`data`](@ref).
"""
properties(img::ImageMeta) = getfield(img, :properties)

Base.hasproperty(img::ImageMeta, k::Symbol) = haskey(properties(img), k)

Base.get(img::ImageMeta, k::Symbol, default) = get(properties(img), k, default)

# So that defaults don't have to be evaluated unless they are needed,
# we also define a @get macro (thanks Toivo Hennington):
struct IMNothing end   # to avoid confusion in the case where dict[key] === nothing
macro get(img, k, default)
    quote
        img, k = $(esc(img)), $(esc(k))
        val = get(properties(img), k, IMNothing())
        return isa(val, IMNothing) ? $(esc(default)) : val
    end
end

ImageAxes.timeaxis(img::ImageMetaAxis) = timeaxis(data(img))
ImageAxes.timedim(img::ImageMetaAxis) = timedim(data(img))
ImageAxes.colordim(img::ImageMetaAxis) = colordim(data(img))

ImageCore.pixelspacing(img::ImageMeta) = pixelspacing(data(img))

"""
    spacedirections(img)

Using ImageMetadata, you can set this property manually. For example, you
could indicate that a photograph was taken with the camera tilted
30-degree relative to vertical using

```
img["spacedirections"] = ((0.866025,-0.5),(0.5,0.866025))
```

If not specified, it will be computed from `pixelspacing(img)`, placing the
spacing along the "diagonal".  If desired, you can set this property in terms of
physical units, and each axis can have distinct units.
"""
ImageCore.spacedirections(img::ImageMeta) = @get img :spacedirections spacedirections(data(img))

ImageCore.sdims(img::ImageMetaAxis) = sdims(data(img))

ImageCore.coords_spatial(img::ImageMetaAxis) = coords_spatial(data(img))

ImageCore.spatialorder(img::ImageMetaAxis) = spatialorder(data(img))

ImageAxes.nimages(img::ImageMetaAxis) = nimages(data(img))

ImageCore.size_spatial(img::ImageMetaAxis) = size_spatial(data(img))

ImageCore.indices_spatial(img::ImageMetaAxis) = indices_spatial(data(img))

ImageCore.assert_timedim_last(img::ImageMetaAxis) = assert_timedim_last(data(img))

#### Permutations over dimensions ####

"""
    permutedims(img, perm, [spatialprops])

When permuting the dimensions of an ImageMeta, you can optionally
specify that certain properties are spatial and they will also be
permuted. `spatialprops` defaults to `spatialproperties(img)`.
"""
permutedims

function permutedims_props!(ret::ImageMeta, ip, spatialprops=spatialproperties(ret))
    if !isempty(spatialprops)
        for prop in spatialprops
            if hasproperty(ret, prop)
                a = getproperty(ret, prop)
                if isa(a, AbstractVector)
                    setproperty!(ret, prop, a[ip])
                elseif isa(a, Tuple)
                    setproperty!(ret, prop, a[ip])
                elseif isa(a, AbstractMatrix) && size(a,1) == size(a,2)
                    setproperty!(ret, prop, a[ip,ip])
                else
                    error("Do not know how to handle property ", prop)
                end
            end
        end
    end
    ret
end

function permutedims(img::ImageMetaAxis, perm)
    p = AxisArrays.permutation(perm, axisnames(data(img)))
    ip = sortperm([p...][[coords_spatial(img)...]])
    permutedims_props!(copyproperties(img, permutedims(data(img), p)), ip)
end
function permutedims(img::ImageMeta, perm)
    ip = sortperm([perm...][[coords_spatial(img)...]])
    permutedims_props!(copyproperties(img, permutedims(data(img), perm)), ip)
end

# Note: `adjoint` does not recurse into ImageMeta properties.
function Base.adjoint(img::ImageMeta{T,2}) where {T<:Real}
    ip = sortperm([2,1][[coords_spatial(img)...]])
    permutedims_props!(copyproperties(img, adjoint(data(img))), ip)
end

function Base.adjoint(img::ImageMeta{T,1}) where T<:Real
    check_empty_spatialproperties(img)
    copyproperties(img, data(img)')
end

"""
    spatialproperties(img)

Return a vector of strings, containing the names of properties that
have been declared "spatial" and hence should be permuted when calling
`permutedims`.  Declare such properties like this:

    img[:spatialproperties] = [:spacedirections]
"""
spatialproperties(img::ImageMeta) = @get img :spatialproperties [:spacedirections]

function check_empty_spatialproperties(img)
    sp = spatialproperties(img)
    for prop in sp
        if hasproperty(img, prop)
            error("spatialproperties must be empty, have $prop")
        end
    end
    nothing
end

#### Low-level utilities ####
function showdictlines(io::IO, dict::AbstractDict, suppress::Set)
    for (k, v) in dict
        if k === :suppress
            continue
        end
        if !in(k, suppress)
            print(io, "\n    ", k, ": ")
            print(IOContext(io, :compact => true), v)
        else
            print(io, "\n    ", k, ": <suppressed>")
        end
    end
end
showdictlines(io::IO, dict::AbstractDict, prop::Symbol) = showdictlines(io, dict, Set([prop]))

# printdictval(io::IO, v) = print(io, v)
# function printdictval(io::IO, v::Vector)
#     for i = 1:length(v)
#         print(io, " ", v[i])
#     end
# end

# converts keyword argument to a dictionary
function kwargs2dict(kwargs)
    d = Dict{Symbol,Any}()
    for (k, v) in kwargs
        d[k] = v
    end
    return d
end

include("operators.jl")
include("deprecations.jl")

end # module

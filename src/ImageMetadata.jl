__precompile__()

module ImageMetadata

# The order here is designed to avoid an ambiguity warning in convert,
# see the top of ImageAxes
using ImageAxes
using ImageCore, Colors, FixedPointNumbers
using ColorVectorSpace   # for overriding math operations with Gray/RGB
using Compat

import Base: +, .+, -, .-, *, .*, /, ./, .^, .<, .>, .==
import Base: +, -, *, /
import Base: permutedims

if VERSION < v"0.6.0-dev.2068"
    const ViewIndex = Base.ViewIndex
else
    const ViewIndex = Union{Base.ViewIndex, Colon}
end

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
`props`), or with `Image(A, prop1=val1, prop2=val2, ...)`.
"""
type ImageMeta{T,N,A<:AbstractArray} <: AbstractArray{T,N}
    data::A
    properties::Dict{String,Any}

    function (::Type{ImageMeta{T,N,A}}){T,N,A}(data::AbstractArray, properties::Dict)
        new{T,N,A}(data, properties)
    end
end
ImageMeta{T,N}(data::AbstractArray{T,N}, props::Dict) = ImageMeta{T,N,typeof(data)}(data,props)
ImageMeta(data::AbstractArray; kwargs...) = ImageMeta(data, kwargs2dict(kwargs))

@compat const ImageMetaArray{T,N,A<:Array} = ImageMeta{T,N,A}
@compat const ImageMetaAxis{T,N,A<:AxisArray} = ImageMeta{T,N,A}

Base.size(A::ImageMeta) = size(A.data)

datatype{T,N,A<:AbstractArray}(::Type{ImageMeta{T,N,A}}) = A

@compat Base.IndexStyle{M<:ImageMeta}(::Type{M}) = IndexStyle(datatype(M))

# getindex and setindex!
for AType in (ImageMeta, ImageMetaAxis)
    @eval begin
        @inline function Base.getindex{T}(img::$AType{T,1}, i::Int)
            @boundscheck checkbounds(img.data, i)
            @inbounds ret = img.data[i]
            ret
        end
        @inline function Base.getindex(img::$AType, i::Int)
            @boundscheck checkbounds(img.data, i)
            @inbounds ret = img.data[i]
            ret
        end
        @inline function Base.getindex{T,N}(img::$AType{T,N}, I::Vararg{Int,N})
            @boundscheck checkbounds(img.data, I...)
            @inbounds ret = img.data[I...]
            ret
        end

        @inline function Base.setindex!{T}(img::$AType{T,1}, val, i::Int)
            @boundscheck checkbounds(img.data, i)
            @inbounds img.data[i] = val
            val
        end
        @inline function Base.setindex!(img::$AType, val, i::Int)
            @boundscheck checkbounds(img.data, i)
            @inbounds img.data[i] = val
            val
        end
        @inline function Base.setindex!{T,N}(img::$AType{T,N}, val, I::Vararg{Int,N})
            @boundscheck checkbounds(img.data, I...)
            @inbounds img.data[I...] = val
            val
        end
    end
end

@inline function Base.getindex(img::ImageMetaAxis, ax::Axis, I...)
    copyproperties(img, img.data[ax, I...])
end
@inline function Base.getindex(img::ImageMetaAxis, i::Union{Integer,AbstractVector,Colon}, I...)
    copyproperties(img, img.data[i, I...])
end

@inline function Base.setindex!(img::ImageMetaAxis, val, ax::Axis, I...)
    setindex!(img.data, val, ax, I...)
end
@inline function Base.setindex!(img::ImageMetaAxis, val, i::Union{Integer,AbstractVector,Colon}, I...)
    setindex!(img.data, val, i, I...)
end

Base.view(img::ImageMeta, ax::Axis, I...) = shareproperties(img, view(img.data, ax, I...))
Base.view{T,N}(img::ImageMeta{T,N}, I::Vararg{ViewIndex,N}) = shareproperties(img, view(img.data, I...))
Base.view(img::ImageMeta, i::ViewIndex) = shareproperties(img, view(img.data, i))
Base.view{N}(img::ImageMeta, I::Vararg{ViewIndex,N}) = shareproperties(img, view(img.data, I...))

Base.getindex(img::ImageMeta, propname::AbstractString) = img.properties[propname]

Base.setindex!(img::ImageMeta, X, propname::AbstractString) = setindex!(img.properties, X, propname)

Base.copy(img::ImageMeta) = ImageMeta(copy(img.data), deepcopy(img.properties))

Base.convert(::Type{ImageMeta}, A::ImageMeta) = A
Base.convert(::Type{ImageMeta}, A::AbstractArray) = ImageMeta(A)
Base.convert{T}(::Type{ImageMeta{T}}, A::ImageMeta{T}) = A
Base.convert{T}(::Type{ImageMeta{T}}, A::ImageMeta) = shareproperties(A, convert(Array{T}, A.data))
Base.convert{T}(::Type{ImageMeta{T}}, A::AbstractArray) = ImageMeta(convert(Array{T}, A))

# copy properties
function Base.copy!(imgdest::ImageMeta, imgsrc::ImageMeta, prop1::AbstractString, props::AbstractString...)
    imgdest[prop1] = deepcopy(imgsrc[prop1])
    for p in props
        imgdest[p] = deepcopy(imgsrc[p])
    end
    imgdest
end

# similar
Base.similar{T}(img::ImageMeta, ::Type{T}, shape::Dims) = ImageMeta(similar(img.data, T, shape), deepcopy(img.properties))
Base.similar{T}(img::ImageMetaAxis, ::Type{T}) = ImageMeta(similar(img.data, T), deepcopy(img.properties))

"""
    copyproperties(img::ImageMeta, data) -> imgnew

Create a new "image," copying the properties dictionary of `img` but
using the data of the AbstractArray `data`. Note that changing the
properties of `imgnew` does not affect the properties of `img`.

See also: [`shareproperties`](@ref).
"""
copyproperties(img::ImageMeta, data::AbstractArray) =
    ImageMeta(data, deepcopy(img.properties))

"""
    shareproperties(img::ImageMeta, data) -> imgnew

Create a new "image," reusing the properties dictionary of `img` but
using the data of the AbstractArray `data`. The two images have
synchronized properties; modifying one also affects the other.

See also: [`copyproperties`](@ref).
"""
shareproperties(img::ImageMeta, data::AbstractArray) = ImageMeta(data, img.properties)

# Delete a property!
Base.delete!(img::ImageMeta, propname::AbstractString) = delete!(img.properties, propname)

# Iteration
# Defer to the array object in case it has special iteration defined
Base.start(img::ImageMeta) = start(data(img))
Base.next(img::ImageMeta, s) = next(data(img), s)
Base.done(img::ImageMeta, s) = done(data(img), s)

# Show
const emptyset = Set()
function showim(io::IO, img::ImageMeta)
    IT = typeof(img)
    print(io, eltype(img).name.name, " ImageMeta with:\n  data: ", summary(img.data), "\n  properties:")
    showdictlines(io, img.properties, get(img, "suppress", emptyset))
end
Base.show(io::IO, img::ImageMeta) = showim(io, img)
Base.show(io::IO, ::MIME"text/plain", img::ImageMeta) = showim(io, img)

Base.reinterpret{T}(::Type{T}, img::ImageMetaArray) = shareproperties(img, reinterpret(T, img.data))
Base.reinterpret{T}(::Type{T}, img::ImageMeta) = error("reinterpret method not defined for $(typeof(img)). Consider a MappedArray instead.")

"""
    data(img::ImageMeta) -> array

Extract the data from `img`, omitting the properties
dictionary. `array` shares storage with `img`, so changes to one
affect the other.

See also: [`properties`](@ref).
"""
ImageCore.data(img::ImageMeta) = img.data   # fixme when deprecation is removed from ImageCore

function ImageCore.permuteddimsview(A::ImageMeta, perm)
    ip = sortperm([perm...][[coords_spatial(A)...]])  # the inverse spatial permutation
    permutedims_props!(copyproperties(A, permuteddimsview(A.data, perm)), ip)
end
ImageCore.channelview(A::ImageMeta) = shareproperties(A, channelview(A.data))
ImageCore.colorview{C<:Colorant}(::Type{C}, A::ImageMeta) = shareproperties(A, colorview(C, A.data))
ImageCore.rawview{T<:Real}(A::ImageMeta{T}) = shareproperties(A, rawview(A.data))
ImageCore.normedview{T<:FixedPoint,S<:Unsigned}(::Type{T}, A::ImageMeta{S}) = shareproperties(A, normedview(T, A.data))

# AxisArrays functions
AxisArrays.axes(img::ImageMetaAxis) = axes(img.data)
AxisArrays.axisdim(img::ImageMetaAxis, ax) = axisdim(img.data, ax)
AxisArrays.axisnames(img::ImageMetaAxis) = axisnames(img.data)
AxisArrays.axisvalues(img::ImageMetaAxis) = axisvalues(img.data)

#### Properties ####

"""
    properties(imgmeta) -> props

Extract the properties dictionary `props` for `imgmeta`. `props`
shares storage with `img`, so changes to one affect the other.

See also: [`data`](@ref).
"""
properties(img::ImageMeta) = img.properties

Base.haskey(img::ImageMeta, k::AbstractString) = haskey(img.properties, k)

Base.get(img::ImageMeta, k::AbstractString, default) = get(img.properties, k, default)

# So that defaults don't have to be evaluated unless they are needed,
# we also define a @get macro (thanks Toivo Hennington):
macro get(img, k, default)
    quote
        img, k = $(esc(img)), $(esc(k))
        index = Base.ht_keyindex(img.properties, k)
        (index > 0) ? img.properties.vals[index] : $(esc(default))
    end
end

ImageAxes.timeaxis(img::ImageMetaAxis) = timeaxis(data(img))
ImageAxes.timedim(img::ImageMetaAxis) = timedim(data(img))
ImageCore.colordim(img::ImageMetaAxis) = colordim(data(img))

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
ImageCore.spacedirections(img::ImageMeta) = @get img "spacedirections" spacedirections(data(img))

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
            if haskey(ret, prop)
                a = ret.properties[prop]
                if isa(a, AbstractVector)
                    ret.properties[prop] = a[ip]
                elseif isa(a, Tuple)
                    ret.properties[prop] = a[ip]
                elseif isa(a, AbstractMatrix) && size(a,1) == size(a,2)
                    ret.properties[prop] = a[ip,ip]
                else
                    error("Do not know how to handle property ", prop)
                end
            end
        end
    end
    ret
end

function permutedims(img::ImageMetaAxis, perm)
    p = AxisArrays.permutation(perm, axisnames(img.data))
    ip = sortperm([p...][[coords_spatial(img)...]])
    permutedims_props!(copyproperties(img, permutedims(img.data, p)), ip)
end
function permutedims(img::ImageMeta, perm)
    ip = sortperm([perm...][[coords_spatial(img)...]])
    permutedims_props!(copyproperties(img, permutedims(img.data, perm)), ip)
end

Base.ctranspose{T<:Real}(img::ImageMeta{T,2}) = permutedims(img, (2,1))
function Base.ctranspose{T<:Real}(img::ImageMeta{T,1})
    check_empty_spatialproperties(img)
    copyproperties(img, img.data')
end

"""
    spatialproperties(img)

Return a vector of strings, containing the names of properties that
have been declared "spatial" and hence should be permuted when calling
`permutedims`.  Declare such properties like this:

    img["spatialproperties"] = ["spacedirections"]
"""
ImageCore.spatialproperties(img::ImageMeta) = @get img "spatialproperties" ["spacedirections"]

function check_empty_spatialproperties(img)
    sp = spatialproperties(img)
    for prop in sp
        if haskey(img, prop)
            error("spatialproperties must be empty, have $prop")
        end
    end
    nothing
end

#### Low-level utilities ####
function showdictlines(io::IO, dict::Dict, suppress::Set)
    for (k, v) in dict
        if k == "suppress"
            continue
        end
        if !in(k, suppress)
            print(io, "\n    ", k, ": ")
            print(IOContext(io, compact=true), v)
        else
            print(io, "\n    ", k, ": <suppressed>")
        end
    end
end
showdictlines(io::IO, dict::Dict, prop::String) = showdictlines(io, dict, Set([prop]))

# printdictval(io::IO, v) = print(io, v)
# function printdictval(io::IO, v::Vector)
#     for i = 1:length(v)
#         print(io, " ", v[i])
#     end
# end

# converts keyword argument to a dictionary
function kwargs2dict(kwargs)
    d = Dict{String,Any}()
    for (k, v) in kwargs
        d[string(k)] = v
    end
    return d
end

include("operators.jl")

end # module

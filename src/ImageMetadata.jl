module ImageMetadata

# The order here is designed to avoid an ambiguity warning in convert,
# see the top of ImageAxes
using ImageAxes
using ImageCore, Colors

export
    # types
    ImageMeta,

    # functions
    copyproperties,
    data,
    getindexim,
    properties,
    shareproperties,
    spatialproperties,
    viewim

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

    function ImageMeta(data::AbstractArray, properties::Dict)
        check_deprecated_properties(data, properties)
        new(data, properties)
    end
end
ImageMeta{T,N}(data::AbstractArray{T,N}, props::Dict) = ImageMeta{T,N,typeof(data)}(data,props)
ImageMeta(data::AbstractArray; kwargs...) = ImageMeta(data, kwargs2dict(kwargs))

typealias ImageMetaAxis{T,N,A<:AxisArray} ImageMeta{T,N,A}

Base.size(A::ImageMeta) = size(A.data)

Base.linearindexing(A::ImageMeta) = Base.linearindexing(A.data)

# getindex and setindex!
@inline function Base.getindex(img::ImageMeta, i::Int)
    @boundscheck checkbounds(img.data, i)
    @inbounds ret = img.data[i]
    ret
end
@inline function Base.getindex{T,N}(img::ImageMeta{T,N}, I::Vararg{Int,N})
    @boundscheck checkbounds(img.data, I...)
    @inbounds ret = img.data[I...]
    ret
end
@inline function Base.setindex!(img::ImageMeta, val, i::Int)
    @boundscheck checkbounds(img.data, i)
    @inbounds img.data[i] = val
    val
end
@inline function Base.setindex!{T,N}(img::ImageMeta{T,N}, val, I::Vararg{Int,N})
    @boundscheck checkbounds(img.data, I...)
    @inbounds img.data[I...] = val
    val
end

Base.getindex(img::ImageMeta, propname::AbstractString) = img.properties[propname]

Base.setindex!(img::ImageMeta, X, propname::AbstractString) = setindex!(img.properties, X, propname)

Base.copy(img::ImageMeta) = ImageMeta(copy(img.data), deepcopy(img.properties))

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
Base.similar{T}(img::ImageMeta, ::Type{T}, shape::Base.NeedsShaping) = ImageMeta(similar(img.data, T, Base.to_shape(shape)), deepcopy(img.properties))
Base.similar{T}(img::ImageMeta, ::Type{T}, shape::Base.DimsOrInds) = ImageMeta(similar(img.data, T, shape), deepcopy(img.properties))

# Create a new "Image" (could be just an Array) copying the properties
# but replacing the data
copyproperties(img::ImageMeta, data::AbstractArray) =
    ImageMeta(data, deepcopy(img.properties))

# Provide new data but reuse (share) the properties
shareproperties(img::ImageMeta, data::AbstractArray) = ImageMeta(data, img.properties)

# Delete a property!
Base.delete!(img::ImageMeta, propname::AbstractString) = delete!(img.properties, propname)

# getindexim and viewim return an ImageMeta. The first copies the
# properties, the second shares them.
getindexim(img::ImageMeta, I...) = copyproperties(img, img.data[I...])
viewim(img::ImageMeta, I...) = shareproperties(img, view(img.data, I...))

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

ImageCore.data(img::ImageMeta) = img.data   # fixme when deprecation is removed from ImageCore

#### Properties ####

properties(img::ImageMeta) = img.properties

Base.haskey(img::ImageMeta, k::AbstractString) = haskey(img.properties, k)

Base.get(img::ImageMeta, k::AbstractString, default) = get(img.properties, k, default)

# So that defaults don't have to be evaluated unless they are needed,
# we also define a @get macro (thanks Toivo Hennington):
macro get(img, k, default)
    quote
        img, k = $(esc(img)), $(esc(k))
        local val
        if !isa(img, ImageMeta)
            val = $(esc(default))
        else
            index = Base.ht_keyindex(img.properties, k)
            val = (index > 0) ? img.properties.vals[index] : $(esc(default))
        end
        val
    end
end

ImageAxes.timedim(img::ImageMetaAxis) = timedim(data(img))

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
ImageCore.spacedirections(img::ImageMeta) = @get img "spacedirections" ImageCore._spacedirections(img)
ImageCore.spacedirections(img::ImageMetaAxis) = @get img "spacedirections" spacedirections(data(img))

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
function Base.permutedims(img::ImageMeta, p::Union{Vector{Int}, Tuple{Int,Vararg{Int}}}, spatialprops = spatialproperties(img))
    if length(p) != ndims(img)
        error("The permutation must have length equal to the number of dimensions")
    end
    if issorted(p) && length(p) == ndims(img)
        return copy(img)
    end
    ip = invperm(p)
    ret = copyproperties(img, permutedims(img.data, p))
    if !isempty(spatialprops)
        ip = sortperm([p...])
        for prop in spatialprops
            a = img.properties[prop]
            if isa(a, AbstractVector)
                ret.properties[prop] = a[ip]
            elseif isa(a, AbstractMatrix) && size(a,1) == size(a,2)
                ret.properties[prop] = a[ip,ip]
            else
                error("Do not know how to handle property ", prop)
            end
        end
    end
    ret
end

Base.permutedims(img::AxisArray, pstr::Union{Vector{Symbol},Tuple{Symbol,Vararg{Symbol}}}) = permutedims(img, permutation(axisnames, pstr))
Base.permutedims(img::ImageMeta, pstr::Union{Vector{Symbol},Tuple{Symbol,Vararg{Symbol}}}, spatialprops::Vector = String[]) = permutedims(img, permutation(axisnames, pstr), spatialprops)

Base.ctranspose{T}(img::ImageMeta{T,2}) = permutedims(img, (2,1))

"""
    spatialproperties(img)

Return a vector of strings, containing the names of properties that
have been declared "spatial" and hence should be permuted when calling
`permutedims`.  Declare such properties like this:

    img["spatialproperties"] = ["spacedirections"]
"""
ImageCore.spatialproperties(img::ImageMeta) = @get img "spatialproperties" ["spacedirections"]

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

include("deprecated.jl")

end # module

module ImagesMeta

using ImagesCore, ImagesAxes

export
    # types
    ImageMeta,

    # functions
    assert_timedim_last,
    assert_xfirst,
    assert_yfirst,
    coords_spatial,
    copyproperties,
    data,
    getindexim,
    height,
    isxfirst,
    isyfirst,
    nimages,
    pixelspacing,
    properties,
    sdims,
    size_spatial,
    shareproperties,
    spacedirections,
    spatialorder,
    spatialproperties,
    timedim,
    width,
    widthheight

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

# Adding or changing a property via setindex!
Base.setindex!(img::AbstractImage, X, propname::AbstractString) = setindex!(img.properties, X, propname)

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
Base.similar{T}(img::ImageMeta, ::Type{T}, shape) = ImageMeta(similar(img.data, T, shape), deepcopy(img.properties))

# Create a new "Image" (could be just an Array) copying the properties
# but replacing the data
copyproperties(img::AbstractArray, data::AbstractArray) = data

copyproperties(img::ImageMeta, data::AbstractArray) =
    ImageMeta(data, deepcopy(img.properties))

# Provide new data but reuse (share) the properties
shareproperties(img::AbstractArray, data::AbstractArray) = data

shareproperties(img::ImageMeta, data::AbstractArray) = ImageMeta(data, img.properties)

# Delete a property!
Base.delete!(img::ImageMeta, propname::AbstractString) = delete!(img.properties, propname)

# getindexim and viewim return an ImageMeta. The first copies the
# properties, the second shares them.
getindexim(img::AbstractImage, I...) = copyproperties(img, img.data[I...])
viewim(img::AbstractImage, I...) = shareproperties(img, view(img.data, I...))

# Iteration
# Defer to the array object in case it has special iteration defined
Base.start(img::ImageMeta) = start(data(img))
Base.next{T,N}(img::ImageMeta{T,N}, s::Tuple{Bool,Base.IteratorsMD.CartesianIndex{N}}) = next(data(img), s)
Base.done{T,N}(img::ImageMeta{T,N}, s::Tuple{Bool,Base.IteratorsMD.CartesianIndex{N}}) = done(data(img), s)
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

data(img::ImageMeta) = img.data

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
        if !isa(img, AbstractImage)
            val = $(esc(default))
        else
            index = Base.ht_keyindex(img.properties, k)
            val = (index > 0) ? img.properties.vals[index] : $(esc(default))
        end
        val
    end
end

"""
    timedim(img) -> d::Int

Return the dimension of the array used for encoding time, or 0 if not
using an axis for this purpose.

Note: if you want to recover information about the time axis, it is
generally much better to use `timeaxis`.
"""
timedim{T,N}(img::AxisArray{T,N}) = _timedim(filter_time_axis(axes(img), ntuple(identity, Val{N})))
timedim(img::ImageMetaAxis) = timedim(data(img))
_timedim(dim::Tuple{Int}) = dim[1]
_timedim(::Tuple{}) = 0

pixelspacing(img::AbstractArray) = ones(ndims(img))
pixelspacing(img::AxisArray) = map(step, axisvalues(img))
pixelspacing(img::ImageMetaAxis) = pixelspacing(data(img))


spacedirections(img::ImageMeta) = @get img "spacedirections" _spacedirections(img)
function _spacedirections(img::AbstractArray)
    ps = pixelspacing(img)
    T = eltype(ps)
    nd = length(ps)
    Vector{T}[(tmp = zeros(T, nd); tmp[i] = ps[i]; tmp) for i = 1:nd]
end

# number of spatial dimensions in the image
@traitfn sdims{AA<:AxisArray; !HasTimeAxis{AA}}(img::AA) = ndims(img)
@traitfn sdims{AA<:AxisArray;  HasTimeAxis{AA}}(img::AA) = ndims(img)-1
sdims(img::ImageMetaAxis) = sdims(data(img))
sdims(img::AbstractArray) = ndims(img)

"""
   coords_spatial(img)

Return a tuple listing the spatial dimensions of `img`.

Note that a better strategy may be to use ImagesAxes and take slices along the time axis.
"""
coords_spatial(img) = collect(1:ndims(img))
coords_spatial{T,N}(img::AxisArray{T,N}) = filter_spatial_axes(axes(img), ntuple(identity, Val{N}))
coords_spatial(img::ImageMetaAxis) = coords_spatial(data(img))

# order of spatial dimensions
spatialorder(img::AxisArray) = filter_spatial_axes(axes(img), axisnames(img))
spatialorder(img::ImageMetaAxis) = spatialorder(data(img))

# number of time slices
nimages(img::AbstractArray) = 1
nimages(img::AxisArray) = _nimages(timeaxis(img))
nimages(img::ImageMetaAxis) = nimages(data(img))
_nimages() = 1
_nimages(ax::Axis) = length(ax)

# size of the spatial grid
size_spatial(img) = size(img)
size_spatial(img::AxisArray) = filter_spatial_axes(axes(img), size(img))
size_spatial(img::ImageMetaAxis) = size_spatial(data(img))

indices_spatial(img) = indices(img)
indices_spatial(img::AxisArray) = filter_spatial_axes(axes(img), indices(img))
indices_spatial(img::ImageMetaAxis) = indices_spatial(data(img))

filter_spatial_axes{N}(axes::NTuple{N,Axis}, items::NTuple{N}) =
    _filter_spatial_axes(axes, items)
@inline @traitfn _filter_spatial_axes{Ax<:Axis;  TimeAxis{Ax}}(axes::Tuple{Ax,Vararg{Any}}, items) =
    _filter_spatial_axes(tail(axes), tail(items))
@inline @traitfn _filter_spatial_axes{Ax<:Axis; !TimeAxis{Ax}}(axes::Tuple{Ax,Vararg{Any}}, items) =
    (items[1], _filter_spatial_axes(tail(axes), tail(items))...)

filter_time_axis{N}(axes::NTuple{N,Axis}, items::NTuple{N}) =
    _filter_time_axis(axes, items)
@inline @traitfn _filter_time_axis{Ax<:Axis; !TimeAxis{Ax}}(axes::Tuple{Ax,Vararg{Any}}, items) =
    _filter_time_axis(tail(axes), tail(items))
@inline @traitfn _filter_time_axis{Ax<:Axis;  TimeAxis{Ax}}(axes::Tuple{Ax,Vararg{Any}}, items) =
    (items[1], _filter_time_axis(tail(axes), tail(items))...)


#### Utilities for writing "simple algorithms" safely ####
# If you don't feel like supporting multiple representations, call these

# Check that the time dimension, if present, is last
assert_timedim_last(img::AbstractArray) = nothing
@traitfn function assert_timedim_last{AA<:AxisArray; HasTimeAxis{AA}}(img::AA)
    ax = axes(img)[end]
    istimeaxis(ax) || error("time dimension is not last")
    nothing
end
@traitfn assert_timedim_last{AA<:AxisArray; !HasTimeAxis{AA}}(img::AA) = nothing
assert_timedim_last(img::ImageMetaAxis) = assert_timedim_last(data(img))

# Spatial storage order
isyfirst(img::AbstractArray) = spatialorder(img)[1] == :y
function assert_yfirst(img)
    if !isyfirst(img)
        error("Image must have y as its first dimension")
    end
end
isxfirst(img::AbstractArray) = spatialorder(img)[1] == :x
function assert_xfirst(img::AbstractArray)
    if !isxfirst(img)
        error("Image must have x as its first dimension")
    end
end

#### Permutations over dimensions ####

# TODO: decide about the default storage order!!
widthheight(img::AbstractArray) = size(img,1), size(img,2)

width(img::AbstractArray) = widthheight(img)[1]
height(img::AbstractArray) = widthheight(img)[2]

# Permute the dimensions of an image, also permuting the relevant properties. If you have non-default properties that are vectors or matrices relative to spatial dimensions, include their names in the list of spatialprops.
function permutedims(img::ImageMeta, p::Union{Vector{Int}, Tuple{Vararg{Int}}}, spatialprops::Vector = spatialproperties(img))
    if length(p) != ndims(img)
        error("The permutation must have length equal to the number of dimensions")
    end
    if issorted(p) && length(p) == ndims(img)
        return copy(img)
    end
    ip = invperm(p)
    ret = copyproperties(img, permutedims(img.data, p))
    sd = timedim(img)
    if sd > 0
        p = setdiff(p, sd)
    end
    if !isempty(spatialprops)
        ip = sortperm(p)
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

permutedims{S<:AbstractString}(img::AbstractArray, pstr::Union{Vector{S}, Tuple{Vararg{S}}}, spatialprops::Vector = spatialproperties(img)) = error("not supported, please switch to ImagesAxes")
permutedims(img::AxisArray, pstr::Union{Vector{S},Tuple{Vararg{S}}}) = permutedims(img, permutation(axisnames, pstr))
permutedims(img::ImageMeta, pstr::Union{Vector{S},Tuple{Vararg{S}}}, spatialprops::Vector = String[]) = permutedims(img, permutation(axisnames, pstr), spatialprops)

# Define the transpose of a 2d image
ctranspose{T}(img::ImageMeta{T,2}) = permutedims(img, (2,1))

# Default list of spatial properties possessed by an image
spatialproperties(img::ImageMeta) = @get img "spatialproperties" String[]

#### Low-level utilities ####
function permutation(to, from)
    n = length(to)
    nf = length(from)
    d = Dict([(from[i], i) for i = 1:length(from)])
    ind = Array(Int, max(n, nf))
    for i = 1:n
        ind[i] = get(d, to[i], 0)
    end
    ind[n+1:nf] = n+1:nf
    ind
end

function showdictlines(io::IO, dict::Dict, suppress::Set)
    for (k, v) in dict
        if k == "suppress"
            continue
        end
        if !in(k, suppress)
            print(io, "\n    ", k, ": ")
            printdictval(io, v)
        else
            print(io, "\n    ", k, ": <suppressed>")
        end
    end
end

printdictval(io::IO, v) = print(io, v)
function printdictval(io::IO, v::Vector)
    for i = 1:length(v)
        print(io, " ", v[i])
    end
end

# converts keyword argument to a dictionary
function kwargs2dict(kwargs)
    d = Dict{ASCIIString,Any}()
    for (k, v) in kwargs
        d[string(k)] = v
    end
    return d
end

end # module

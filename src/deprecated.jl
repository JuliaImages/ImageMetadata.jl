using IndirectArrays

using Base: depwarn
import Base: sqrt, atan2, hypot, real, imag, abs

#### Types and constructors ####

Base.@deprecate_binding Image ImageMeta
Base.@deprecate_binding AbstractImage ImageMeta
Base.@deprecate_binding AbstractImageDirect ImageMeta
typealias ImageMetaIndirect{T,N,A<:IndirectArray} ImageMeta{T,N,A}
Base.@deprecate_binding AbstractImageIndexed ImageMetaIndirect

@deprecate ImageCmap(data, cmap; kwargs...)  ImageMeta(IndirectArray(data, cmap); kwargs...)
@deprecate ImageCmap(data, cmap, properties) ImageMeta(IndirectArray(data, cmap), properties)

Base.@deprecate_binding sliceim viewim

function subim(img::Union{AxisArray,ImageMeta}, args...)
    newargs = _subim_indexes(args)
    newargstr = join(map(string, newargs), ", ")
    Base.depwarn("subim is deprecated, call viewim(img, $newargstr) instead", :subim)
    viewim(img, newargs...)
end
export subim

# This is not type-stable, but since it's used for a deprecation this is fine
function _subim_indexes(args)
    n = length(args)
    newargs = Array{Any}(n)
    haveextended = false
    for i = n:-1:1
        a = args[i]
        if isa(a, Real)
            newargs[i] = haveextended ? (a:a) : a
        else
            newargs[i] = a
            haveextended |= isa(a, AbstractArray)
        end
    end
    newargs
end

#### Indexing ####

using ImageAxes: getaxes

function Base.view(img::ImageMetaAxis, dimname::AbstractString, ind::Base.ViewIndex, args...)
    axs = getaxes(dimname, ind, args...)
    Base.depwarn("indexing with strings is deprecated, use view(img, $(axs...)) instead", :view!)
    view(img.data, axs...)
end

function viewim(img::ImageMetaAxis, dimname::AbstractString, ind::Base.ViewIndex, args...)
    axs = getaxes(dimname, ind, args...)
    Base.depwarn("indexing with strings is deprecated, use view(img, $(axs...)) instead", :view!)
    shareproperties(img, view(img.data, axs...))
end

@deprecate copyproperties(img::AbstractArray, data::AbstractArray) data
@deprecate shareproperties(img::AbstractArray, data::AbstractArray) data

@deprecate getindexim(img::AbstractArray, I...) img[I...]
@deprecate viewim(img::AbstractArray, I...) view(img, I...)

#### Properties ####

properties_depwarn() = depwarn("""
properties is deprecated for arrays that are not ImageMeta; perhaps use keyword arguments
instead, structuring your code like this:
    function myfunction(img::AbstractArray; pixelspacing=(1,1), ...)
        # code
    end
    myfunction(img::ImageMeta) = myfunction(data(img), pixelspacing=img["pixelspacing"])
""", :properties)
function properties(A::AbstractArray)
    properties_depwarn()
    Dict{String,Any}()
end

function check_deprecated_properties(data, properties)
    for key in ("colorspace", "colordim")
        if haskey(properties, key)
            error("\"$key\" property is ignored, now color is encoded only by the element type.\nSee `colorview` to represent a numeric array as a color array.")
        end
    end
    if haskey(properties, "limits")
        local z, o
        try
            T = eltype(data)
            z = zero(T)
            o = mapc(x->one(x), z)
        catch
            z, o = 0, 1
        end
        error("\"limits\" property is ignored, limits are always ($z,$o)")
    end
    if haskey(properties, "timedim")
        error("\"timedim\" property is ignored, please use the AxisArrays package and A = AxisArray(data, ..., :time, ...)")
    end
    if haskey(properties, "pixelspacing")
        ps = properties["pixelspacing"]
        error("\"pixelspacing\" property is ignored, please use the AxisArrays package and A = AxisArray(data, (axisnames...), ($(join(ps, ','))))")
    end
    if haskey(properties, "spatialorder")
        so = properties["spatialorder"]
        try
            so = map(s->":"*s, so)
        catch
        end
        error("\"spatialorder\" property is ignored, please use the AxisArrays package and A = AxisArray(data, $(join(so, ", ")))")
    end
    nothing
end

import Base: haskey, get
@deprecate haskey(a::AbstractArray, k::AbstractString) false

@deprecate get(img::AbstractArray, k::AbstractString, default) default

### traits ###

import ImageCore.isdirect
@deprecate isdirect(img::IndirectArray) false
@deprecate isdirect(img::ImageMeta) isdirect(data(img))

# function pixelspacing{T,N}(img::ImageMeta{T,N})
#     depwarn("pixelspacing(::ImageMeta) is deprecated, please switch to ImageAxes", :pixelspacing)
#     _pixelspacing(img)
# end
# _pixelspacing{T}(img::ImageMeta{T}) = @get img "pixelspacing" __pixelspacing(img)
# function __pixelspacing(img::ImageMeta)
#     if haskey(img, "spacedirections")
#         sd = img["spacedirections"]
#         return ___pixelspacing(sd)
#     end
#     ones(sdims(img))
# end

# @noinline ___pixelspacing(sd) = [maximum(abs(sd[i])) for i = 1:length(sd)]

# This is mostly for user information---in code it's generally better
# to use spatialorder, colordim, and timedim directly
import ImageAxes.storageorder
@deprecate storageorder(A::AxisArray) axisnames(A)


# Elementwise functions
@deprecate sqrt(img::ImageMeta) shareproperties(img, sqrt.(data(img)))
@deprecate atan2(img1::ImageMeta, img2::ImageMeta) shareproperties(img1, atan2.(data(img1),data(img2)))
@deprecate hypot(img1::ImageMeta, img2::ImageMeta) shareproperties(img1, hypot.(data(img1),data(img2)))
@deprecate real(img::ImageMeta) shareproperties(img,real.(data(img)))
@deprecate imag(img::ImageMeta) shareproperties(img,imag.(data(img)))
@deprecate abs(img::ImageMeta) shareproperties(img,abs.(data(img)))

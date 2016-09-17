using IndirectArrays

import Base: sqrt, atan2, hypot, real, imag, abs

#### Types and constructors ####

Base.@deprecate_binding Image ImageMeta
Base.@deprecate_binding AbstractImage ImageMeta
Base.@deprecate_binding AbstractImageDirect ImageMeta
typealias ImageMetaIndirect{T,N,A<:IndirectArray} ImageMeta{T,N,A}
Base.@deprecate_binding AbstractImageIndexed ImageMetaIndirect

@deprecate ImageCmap(data, cmap; kwargs...)  ImageMeta(IndirectArray(data, cmap); kwargs...)
@deprecate ImageCmap(data, cmap, properties) ImageMeta(IndirectArray(data, cmap), properties)

Base.@deprecate_binding subim viewim
Base.@deprecate_binding sliceim viewim

#### Indexing ####

function getaxes(dimname::AbstractString, ind, nameind...)
    ax1 = Axis{Symbol(dimname)}(ind)
    axs = []
    for i = 1:2:length(nameind)
        push!(axs, Axis{Symbol(nameind[i])}(nameind[i+1]))
    end
    (ax1, axs...)
end

throw_axisarray(dimname, ind) = error("for named dimensions, please index as img[Axis{:$dimname}($ind), ...]")
throw_axisarray(dimname, ind, fsym) = error("for named dimensions, please use $fsym(img, [Axis{:$dimname}($ind), ...)")

function Base.getindex(img::AxisArray, dimname::AbstractString, ind::Base.ViewIndex, nameind...)
    axs = getaxes(dimname, ind, nameind...)
    Base.depwarn("indexing with strings is deprecated, use img[$(axs...)] instead", :setindex!)
    img[axs...]
end

function Base.setindex!(img::AxisArray, X, dimname::AbstractString, ind::Base.ViewIndex, nameind...)
    axs = getaxes(dimname, ind, nameind...)
    Base.depwarn("indexing with strings is deprecated, use img[$(axs...)] instead", :setindex!)
    setindex!(img, X, axs...)
end
Base.setindex!(img::AbstractArray, X, dimname::AbstractString, ind::Base.ViewIndex, nameind...) = error("for named dimensions, please switch to ImageAxes")

Base.view(img::ImageMeta, dimname::AbstractString, ind::Base.ViewIndex, args...) = error("for named dimensions, please switch to ImageAxes")
function Base.view(img::Union{AxisArray, ImageMetaAxis}, dimname::AbstractString, ind::Base.ViewIndex, args...)
    axs = getaxes(dimname, ind, args...)
    Base.depwarn("indexing with strings is deprecated, use view(img, $(axs...)) instead", :view!)
    view(img, axs...)
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
    Dict("colorspace" => colorspace(A),
         "colordim" => colordim(A),
         "timedim" => timedim(A),
         "pixelspacing" => pixelspacing(A),
         "spatialorder" => spatialorder(A))
end
function properties{C<:Colorant}(A::AbstractArray{C})
    properties_depwarn()
    Dict("timedim" => timedim(A),
         "pixelspacing" => pixelspacing(A),
         "spatialorder" => spatialorder(A))
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
            z, o = zero(eltype(data)), one(eltype(data))
        catch
            z, o = 0, 1
        end
        error("\"limits\" property is ignored, limits are always ($z,$o)")
    end
    if haskey(properties, "timedim")
        error("\"timedim\" property is ignored, please switch to ImageAxes and use A = AxisArray(data, ..., :time, ...)")
    end
    if haskey(properties, "pixelspacing")
        error("\"pixelspacing\" property is ignored, please switch to ImageAxes and use A = AxisArray(data, ..., Axis{:y}(start:[step:]stop), ...)")
    end
    if haskey(properties, "spatialorder")
        so = properties["spatialorder"]
        try
            so = map(s->":"*s, so)
        catch
        end
        error("\"spatialorder\" property is ignored, please switch to ImageAxes and use A = AxisArray(data, $(join(so, ", ")))")
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

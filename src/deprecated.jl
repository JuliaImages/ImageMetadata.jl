using IndirectArrays

#### Types and constructors ####

# Perhaps move these to Images.jl?
Base.@deprecate_binding Image ImageMeta
Base.@deprecate_binding AbstractImage ImageMeta
Base.@deprecate_binding AbstractImageDirect ImageMeta
typealias ImageMetaIndirect{T,N,A<:IndirectArray} ImageMeta{T,N,A}
Base.@deprecate_binding AbstractImageIndexed ImageMetaIndirect

@deprecate ImageCmap(data, cmap) IndirectArray(data, cmap)
@deprecate ImageCmap(data, cmap, properties) ImageMeta(IndirectArray(data, cmap), properties)

#### Indexing ####

Base.setindex!(img::AxisArray, X, dimname::AbstractString, ind::Base.ViewIndex, nameind...) = error("for named dimensions, please switch to ImagesAxes")
Base.setindex!(img::AbstractArray, X, dimname::AbstractString, ind::Base.ViewIndex, nameind...) = error("for named dimensions, please switch to ImagesAxes")

Base.view(img::ImageMeta, dimname::AbstractString, ind::Base.ViewIndex, args...) = error("for named dimensions, please switch to ImagesAxes")

@deprecate copyproperties(img::AbstractArray, data::AbstractArray) data
@deprecate shareproperties(img::AbstractArray, data::AbstractArray) data

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
        error("\"timedim\" property is ignored, please switch to ImagesAxes and use A = AxisArray(data, ..., :time, ...)")
    end
    if haskey(properties, "pixelspacing")
        error("\"pixelspacing\" property is ignored, please switch to ImagesAxes and use A = AxisArray(data, ..., Axis{:y}(start:[step:]stop), ...)")
    end
    if haskey(properties, "spatialorder")
        so = properties["spatialorder"]
        try
            so = map(s->":"*s, so)
        catch
        end
        error("\"spatialorder\" property is ignored, please switch to ImagesAxes and use A = AxisArray(data, $(join(so, ", ")))")
    end
    nothing
end

import Base: haskey, get
@deprecate haskey(a::AbstractArray, k::AbstractString) false

@deprecate get(img::AbstractArray, k::AbstractString, default) default

### traits ###

import ImagesCore.isdirect
@deprecate isdirect(img::IndirectArray) false
@deprecate isdirect(img::ImageMeta) isdirect(data(img))

# function pixelspacing{T,N}(img::ImageMeta{T,N})
#     depwarn("pixelspacing(::ImageMeta) is deprecated, please switch to ImagesAxes", :pixelspacing)
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
import ImagesAxes.storageorder
@deprecate storageorder(A::AxisArray) axisnames(A)

#### Permutations over dimensions ####

import Base.permutedims
@deprecate permutedims{S<:AbstractString}(img::AxisArray, pstr::Union{Vector{S},Tuple{S,Vararg{S}}}, spatialprops::Vector) permutedims(img, pstr)

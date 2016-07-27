using IndirectArrays

#### Types and constructors ####

# Probably move these to Images.jl:
Base.@deprecate_binding Image ImageMeta

@deprecate ImageCmap(data, cmap) IndirectArray(data, cmap)
@deprecate ImageCmap(data, cmap, properties) ImageMeta(IndirectArray(data, cmap), properties)

#### Indexing ####

Base.setindex!(img::AbstractImage, X, dimname::AbstractString, ind::Base.ViewIndex, nameind...) = error("for named dimensions, please switch to ImagesAxes")

Base.view(img::AbstractImage, dimname::AbstractString, ind::Base.ViewIndex, args...) = error("for named dimensions, please switch to ImagesAxes")

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

import Base: haskey, get
@deprecate haskey(a::AbstractArray, k::AbstractString) false

@deprecate get(img::AbstractArray, k::AbstractString, default) default

### traits ###

@deprecate isdirect(img::IndirectArray) false

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
@deprecate storageorder(A::AxisArray) axisnames(A)

#### Permutations over dimensions ####

@deprecate permutedims{S<:AbstractString}(img::AxisArray, pstr::Union{Vector{S},Tuple{Vararg{S}}}, spatialprops::Vector) permutedims(img, pstr)

using IndirectArrays

#### Types and constructors ####

@deprecate_binding AbstractImage AbstractArray
@deprecate_binding AbstractImageDirect AbstractArray
@deprecate_binding AbstractImageIndexed AbstractArray

@deprecate_binding Image ImageMeta

@deprecate ImageCmap(data, cmap) IndirectArray(data, cmap)
@deprecate ImageCmap(data, cmap, properties) ImageMeta(IndirectArray(data, cmap), properties)

# Convenience constructors
@deprecate grayim(A) ColorView{Gray}(A)
grayim(A::AbstractArray{UInt8,2})  = grayim(reinterpret(UFixed8, A))
grayim(A::AbstractArray{UInt16,2}) = grayim(reinterpret(UFixed16, A))
grayim(A::AbstractArray{UInt8,3})  = grayim(reinterpret(UFixed8, A))
grayim(A::AbstractArray{UInt16,3}) = grayim(reinterpret(UFixed16, A))

export colorim
function colorim{T}(A::AbstractArray{T,3})
    if size(A, 1) == 4 || size(A, 3) == 4
        error("The array looks like a 4-channel color image. Please specify the colorspace explicitly (e.g. \"ARGB\" or \"RGBA\".)")
    end

    colorim(A, "RGB")
end
function colorim{T<:Fractional}(A::AbstractArray{T,3}, colorspace)
    Base.depwarn("colorim(A, colorspace) is deprecated, use ColorView{C}(A) instead, possibly in conjunction with permutedview", :colorim)
    CT = getcolortype(colorspace, eltype(A))
    if 3 <= size(A, 1) <= 4 && 3 <= size(A, 3) <= 4
        error("Both first and last dimensions are of size 3 or 4; impossible to guess which is for color. Use the Image constructor directly.")
    elseif 3 <= size(A, 1) <= 4  # Image as returned by imread for regular 2D RGB images
        ColorView{CT}(A)
    elseif 3 <= size(A, 3) <= 4  # "Matlab"-style image, as returned by convert(Array, im).
        ColorView{CT}(permutedview(A))
    else
        error("Neither the first nor the last dimension is of size 3. This doesn't look like an RGB image.")
    end
end
colorim(A::AbstractArray{UInt8,3},  colorspace) = colorim(reinterpret(UFixed8, A), colorspace)
colorim(A::AbstractArray{UInt16,3}, colorspace) = colorim(reinterpret(UFixed16, A), colorspace)

Base.setindex!(img::AbstractImage, X, dimname::AbstractString, ind::Base.ViewIndex, nameind...) = error("for named dimensions, please switch to ImagesAxes")

Base.view(img::AbstractImage, dimname::AbstractString, ind::Base.ViewIndex, args...) = error("for named dimensions, please switch to ImagesAxes")

function data(img::AbstractArray)
    Base.depwarn("""
data(A) is deprecated for arrays that are not an ImageMeta; structure your code like this:
    function myfunction(img::AbstractArray, args...)
        # "real" algorithm goes here
    end
    myfunction(img::ImageMeta, args...) = myfunction(data(img), args...)
""", :data)
    img
end

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

@deprecate Base.haskey(a::AbstractArray, k::AbstractString) = false

@deprecate Base.get(img::AbstractArray, k::AbstractString, default) = default

# Using plain arrays, we have to make all sorts of guesses about
# colorspace and storage order. This can be a big problem for
# three-dimensional images, image sequences, cameras with more than
# 16-bits, etc. In such cases use an AbstractImage type.

# Here are the two most important assumptions (see also colorspace below):
defaultarraycolordim = 3
# defaults for plain arrays ("vertical-major")
const yx = ["y", "x"]
# order used in Cairo & most image file formats (with color as the very first dimension)
const xy = ["x", "y"]
function spatialorder(img::AbstractArray)
    depwarn("spatialorder is deprecated, please switch to ImagesAxes instead", :spatialorder)
    _spatialorder(img)
end
_spatialorder(::Type{Matrix}) = yx
_spatialorder(img::AbstractArray) = (sdims(img) == 2) ? spatialorder(Matrix) : error("cannot guess spatial order for images with ", sdims(img), " spatial dimensions")

@deprecate isdirect(img::AbstractArray) true
@deprecate isdirect(img::IndirectArray) false

export colorspace
function colorspace(img)
    depwarn("""
colorspace(img) is deprecated, use eltype(img) instead, possibly in conjunction
with colorview(img)""", :colorspace)
    _colorspace(img)
end
_colorspace{C<:Colorant}(img::AbstractVector{C}) = ColorTypes.colorant_string(C)
_colorspace{C<:Colorant}(img::AbstractMatrix{C}) = ColorTypes.colorant_string(C)
_colorspace{C<:Colorant}(img::AbstractArray{C,3}) = ColorTypes.colorant_string(C)
_colorspace{C<:Colorant}(img::AbstractImage{C}) = ColorTypes.colorant_string(C)
_colorspace(img::AbstractVector{Bool}) = "Binary"
_colorspace(img::AbstractMatrix{Bool}) = "Binary"
_colorspace(img::AbstractArray{Bool}) = "Binary"
_colorspace(img::AbstractArray{Bool,3}) = "Binary"
_colorspace(img::AbstractMatrix{UInt32}) = "RGB24"
_colorspace(img::AbstractVector) = "Gray"
_colorspace(img::AbstractMatrix) = "Gray"
_colorspace{T}(img::AbstractArray{T,3}) = (size(img, defaultarraycolordim) == 3) ? "RGB" : error("Cannot infer colorspace of Array, use an AbstractImage type")
_colorspace(img::ImageMeta{Bool}) = "Binary"
csinfer{C<:Colorant}(::Type{C}) = ColorTypes.colorant_string(C)
csinfer(C) = "Unknown"
_colorspace(img::AbstractImage) = get(img.properties, "_colorspace", "Unknown")

colorspacedict = Dict{ASCIIString,Any}()
for ACV in (Color, AbstractRGB)
    for CV in subtypes(ACV)
        (length(CV.parameters) == 1 && !(CV.abstract)) || continue
        str = string(CV.name.name)
        colorspacedict[str] = CV
    end
end
function getcolortype{T}(str::ASCIIString, ::Type{T})
    if haskey(colorspacedict, str)
        CV = colorspacedict[str]
        return CV{T}
    else
        if endswith(str, "A")
            CV = colorspacedict[str[1:end-1]]
            return coloralpha(CV){T}
        elseif startswith(str, "A")
            CV = colorspacedict[str[2:end]]
            return alphacolor(CV){T}
        else
            error("colorspace $str not recognized")
        end
    end
end

export colordim
function colordim(img)
    depwarn("colordim(img) is deprecated, use colorview(img) to represent as a color image", :colordim)
    _colordim(img)
end
_colordim{C<:Colorant}(img::AbstractVector{C}) = 0
_colordim{C<:Colorant}(img::AbstractMatrix{C}) = 0
_colordim{C<:Colorant}(img::AbstractArray{C,3}) = 0
_colordim{C<:Colorant}(img::ImageMeta{C}) = 0
_colordim(img::AbstractVector) = 0
_colordim(img::AbstractMatrix) = 0
_colordim{T}(img::ImageMeta{T,3}) = get(img, "colordim", 0)::Int

export timedim
function timedim(img)
    depwarn("timedim(img) is deprecated, please switch to ImagesAxes and use `timeaxis` instead", :timedim)
    return 0
end

export limits
function limits(img)
    depwarn("limits(img) is deprecated, use (zero(T),one(T)) where T is the eltype", :limits)
    _limits(img)
end
_limits(img::AbstractArray{Bool}) = 0,1
_limits{T<:AbstractFloat}(img::AbstractArray{T}) = zero(T), one(T)
_limits(img::AbstractImage{Bool}) = 0,1
_limits{T}(img::AbstractImageDirect{T}) = get(img, "limits", (zero(T), one(T)))
_limits(img::AbstractImageIndexed) = @get img "limits" (minimum(img.cmap), maximum(img.cmap))

function pixelspacing{T,N}(img::ImageMeta{T,N})
    depwarn("pixelspacing(::ImageMeta) is deprecated, please switch to ImagesAxes", :pixelspacing)
    _pixelspacing(img)
end
_pixelspacing{T}(img::ImageMeta{T}) = @get img "pixelspacing" __pixelspacing(img)
function __pixelspacing(img::ImageMeta)
    if haskey(img, "spacedirections")
        sd = img["spacedirections"]
        return ___pixelspacing(sd)
    end
    ones(sdims(img))
end

@noinline ___pixelspacing(sd) = [maximum(abs(sd[i])) for i = 1:length(sd)]

# This is mostly for user information---in code it's generally better
# to use spatialorder, colordim, and timedim directly
@deprecate storageorder(A::AxisArray) axisnames(A)
function storageorder(img::AbstractArray)
    depwarn("storageorder is deprecated, please switch to ImagesAxes and use `axisnames`", :storageorder)
    so = Array(String, ndims(img))
    so[coords_spatial(img)] = spatialorder(img)
    td = timedim(img)
    if td != 0
        so[td] = "t"
    end
    so
end

# number of array elements used for each pixel/voxel
@deprecate ncolorelem{C<:Colorant}(img::AbstractArray{C}) length(C)
function ncolorelem(img)
    depwarn("ncolorelem is deprecated, please encode as a color array (possibly with `colorview`) and use `length(eltype(img))`.\nNumeric arrays are assumed to be grayscale and will return 1.")
    1
end

#### Utilities for writing "simple algorithms" safely ####
# If you don't feel like supporting multiple representations, call these

# Two-dimensional images
export assert2d
function assert2d(img::AbstractArray)
    depwarn("assert2d is deprecated, write your algorithm as `myfunc{T}(img::AbstractArray{T,2}) instead", :assert2d)
    if ndims(img) != 2
        error("Only two-dimensional images are supported")
    end
end

# "Scalar color", either grayscale, RGB24, or an immutable type
export assert_scalar_color
function assert_scalar_color(img::AbstractArray)
    depwarn("assert_scalar_color is deprecated and can be removed", :assert_scalar_color)
    nothing
end


#### Permutations over dimensions ####

function spatialproperties(img::AbstractArray)
    depwarn("spatialproperties is deprecated for any arrays other than ImageMeta", :spatialproperties)
    String[]
end

@deprecate spatialpermutation permutation

# width and height, translating "x" and "y" spatialorder into horizontal and vertical, respectively
# Permute the dimensions of an image, also permuting the relevant properties. If you have non-default properties that are vectors or matrices relative to spatial dimensions, include their names in the list of spatialprops.
@deprecate permutedims{S<:AbstractString}(img::AbstractArray, pstr::Union{Vector{S}, Tuple{Vararg{S}}}, spatialprops::Vector = spatialproperties(img)) permutedims(img, map(Symbol, pstr), spatialprops)
@deprecate permutedims(img::AxisArray, pstr::Union{Vector{S},Tuple{Vararg{S}}}, spatialprops::Vector) permutedims(img, pstr)

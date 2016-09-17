(+)(img::ImageMeta{Bool}, n::Bool) = img .+ n
(+)(n::Bool, img::ImageMeta{Bool}) = n .+ img
(+)(img::ImageMeta, n::Number) = img .+ n
(+)(img::ImageMeta, n::AbstractRGB) = img .+ n
(+)(n::Number, img::ImageMeta) = n .+ img
(+)(n::AbstractRGB, img::ImageMeta) = n .+ img
(.+)(img::ImageMeta, n::Number) = shareproperties(img, data(img).+n)
(.+)(n::Number, img::ImageMeta) = shareproperties(img, data(img).+n)
if isdefined(:UniformScaling)
    (+){Timg,TA<:Number}(img::ImageMeta{Timg,2}, A::UniformScaling{TA}) = shareproperties(img, data(img)+A)
    (-){Timg,TA<:Number}(img::ImageMeta{Timg,2}, A::UniformScaling{TA}) = shareproperties(img, data(img)-A)
end
(+)(img::ImageMeta, A::BitArray) = shareproperties(img, data(img)+A)
(+)(img::ImageMeta, A::ImageMeta) = shareproperties(img, data(img)+data(A))
(+)(img::ImageMeta, A::AbstractArray) = shareproperties(img, data(img)+data(A))
(+){S,T}(A::Range{S}, img::ImageMeta{T}) = shareproperties(img, data(A)+data(img))
(+)(A::AbstractArray, img::ImageMeta) = shareproperties(img, data(A)+data(img))
(.+)(img::ImageMeta, A::BitArray) = shareproperties(img, data(img).+A)
(.+)(img::ImageMeta, A::AbstractArray) = shareproperties(img, data(img).+data(A))
(-)(img::ImageMeta{Bool}, n::Bool) = img .- n
(-)(img::ImageMeta, n::Number) = img .- n
(-)(img::ImageMeta, n::AbstractRGB) = img .- n
(.-)(img::ImageMeta, n::Number) = shareproperties(img, data(img).-n)
(-)(n::Bool, img::ImageMeta{Bool}) = n .- img
(-)(n::Number, img::ImageMeta) = n .- img
(-)(n::AbstractRGB, img::ImageMeta) = n .- img
(.-)(n::Number, img::ImageMeta) = shareproperties(img, n.-data(img))
(-)(img::ImageMeta, A::BitArray) = shareproperties(img, data(img)-A)
(-){T}(img::ImageMeta{T,2}, A::Diagonal) = shareproperties(img, data(img)-A) # fixes an ambiguity warning
(-)(img::ImageMeta, A::Range) = shareproperties(img, data(img)-A)
(-)(img::ImageMeta, A::ImageMeta) = shareproperties(img, data(img)-data(A))
(-)(img::ImageMeta, A::AbstractArray) = shareproperties(img, data(img)-data(A))
(-){S,T}(A::Range{S}, img::ImageMeta{T}) = shareproperties(img, data(A)-data(img))
(-)(A::AbstractArray, img::ImageMeta) = shareproperties(img, data(A)-data(img))
(-)(img::ImageMeta) = shareproperties(img, -data(img))
(.-)(img::ImageMeta, A::BitArray) = shareproperties(img, data(img).-A)
(.-)(img::ImageMeta, A::AbstractArray) = shareproperties(img, data(img).-data(A))
(*)(img::ImageMeta, n::Number) = (.*)(img, n)
(*)(n::Number, img::ImageMeta) = (.*)(n, img)
(.*)(img::ImageMeta, n::Number) = shareproperties(img, data(img).*n)
(.*)(n::Number, img::ImageMeta) = shareproperties(img, data(img).*n)
(/)(img::ImageMeta, n::Number) = shareproperties(img, data(img)/n)
(.*)(img1::ImageMeta, img2::ImageMeta) = shareproperties(img1, data(img1).*data(img2))
(.*)(img::ImageMeta, A::BitArray) = shareproperties(img, data(img).*A)
(.*)(A::BitArray, img::ImageMeta) = shareproperties(img, data(img).*A)
(.*)(img::ImageMeta{Bool}, A::BitArray) = shareproperties(img, data(img).*A)
(.*)(A::BitArray, img::ImageMeta{Bool}) = shareproperties(img, data(img).*A)
(.*)(img::ImageMeta, A::AbstractArray) = shareproperties(img, data(img).*A)
(.*)(A::AbstractArray, img::ImageMeta) = shareproperties(img, data(img).*A)
(./)(img::ImageMeta, A::BitArray) = shareproperties(img, data(img)./A)  # needed to avoid ambiguity warning
(./)(img1::ImageMeta, img2::ImageMeta) = shareproperties(img1, data(img1)./data(img2))
(./)(img::ImageMeta, A::AbstractArray) = shareproperties(img, data(img)./A)
(.^)(img::ImageMeta, p::Number) = shareproperties(img, data(img).^p)

# Logical operations
(.<)(img::ImageMeta, n::Number) = data(img) .< n
(.>)(img::ImageMeta, n::Number) = data(img) .> n
(.<)(img::ImageMeta{Bool}, A::AbstractArray{Bool}) = data(img) .< A
(.<)(img::ImageMeta, A::AbstractArray) = data(img) .< A
(.>)(img::ImageMeta, A::AbstractArray) = data(img) .> A
(.==)(img::ImageMeta, n::Number) = data(img) .== n
(.==)(img::ImageMeta{Bool}, A::AbstractArray{Bool}) = data(img) .== A
(.==)(img::ImageMeta, A::AbstractArray) = data(img) .== A

fft(img::ImageMeta) = shareproperties(img, fft(data(img)))
function fft(img::ImageMeta, region, args...)
    F = fft(data(img), region, args...)
    props = copy(properties(img))
    props["region"] = region
    ImageMeta(F, props)
end

function ifft(img::ImageMeta)
    region = get(img, "region", 1:ndims(img))
    A = ifft(data(img), region)
    props = copy(properties(img))
    haskey(props, "region") && delete!(props, "region")
    ImageMeta(A, props)
end
ifft(img::ImageMeta, region, args...) = ifft(data(img), region, args...)

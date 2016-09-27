using ColorVectorSpace: AbstractGray, TransparentGray, TransparentRGB

(-)(img::ImageMeta) = shareproperties(img, -data(img))

for op in (:+, :.+, :-, :.-, :*, :.*)
    @eval begin
        ($op)(img::ImageMeta{Bool}, n::Bool) = shareproperties(img, ($op)(data(img), n))
        ($op)(n::Bool, img::ImageMeta{Bool}) = shareproperties(img, ($op)(n, data(img)))
        ($op)(img::ImageMeta, n::Number) = shareproperties(img, ($op)(data(img), n))
        ($op)(n::Number, img::ImageMeta) = shareproperties(img, ($op)(n, data(img)))
    end
    if op !== :*
        @eval begin
            ($op)(img::ImageMeta, B::BitArray) = shareproperties(img, ($op)(data(img), B))
            ($op)(B::BitArray, img::ImageMeta) = shareproperties(img, ($op)(B, data(img)))
            ($op)(img::ImageMeta, B::ImageMeta) = ambigop(Symbol($op))
            ($op)(img::ImageMeta, B::AbstractArray) = shareproperties(img, ($op)(data(img), B))
            ($op)(B::AbstractArray, img::ImageMeta) = shareproperties(img, ($op)(B, data(img)))
        end
    end
end

for op in (:+, :.+, :-, :.-)
    for CV in (:AbstractGray, :TransparentGray, :AbstractRGB, :TransparentRGB)
        @eval begin
            ($op){CV<:$CV}(img::ImageMeta{CV}, n::$CV) =
                shareproperties(img, ($op)(data(img), n))
            ($op){CV<:$CV}(n::$CV, img::ImageMeta{CV}) =
                shareproperties(img, ($op)(n, data(img)))
        end
    end
end

(/)(img::ImageMeta, n::Number) = shareproperties(img, data(img)/n)
(./)(img::ImageMeta, A::BitArray) = shareproperties(img, data(img)./A)  # needed to avoid ambiguity warning
(./)(img1::ImageMeta, img2::ImageMeta) = ambigop(:./)
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

ambigop(s::Symbol) = error("$s with two ImageMeta arrays: dictionary choice is ambiguous")

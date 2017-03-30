using ColorVectorSpace: AbstractGray, TransparentGray, TransparentRGB

const dotops = VERSION < v"0.6.0-dev.1839"

if isdefined(Base.Broadcast, :containertype)
    # Specialize the low-level broadcasting machinery for ImageMeta.
    # This make it work for all operators, rather than specializing
    # each operator individually.
    using Base.Broadcast: _broadcast_eltype, broadcast_indices
    Base.Broadcast._containertype(::Type{<:ImageMeta}) = ImageMeta
    Base.Broadcast.promote_containertype(::Type{ImageMeta}, ::Type{ImageMeta}) = ImageMeta
    Base.Broadcast.promote_containertype(::Type{Array}, ::Type{ImageMeta}) = ImageMeta
    Base.Broadcast.promote_containertype(::Type{ImageMeta}, ::Type{Array}) = ImageMeta
    Base.Broadcast.promote_containertype(::Type{ImageMeta}, ct) = ImageMeta
    Base.Broadcast.promote_containertype(ct, ::Type{ImageMeta}) = ImageMeta
    Base.Broadcast.broadcast_indices(::Type{ImageMeta}, A) = indices(A)
    function Base.Broadcast.broadcast_c(f, ::Type{ImageMeta}, As...)
        T = _broadcast_eltype(f, As...)
        shape = broadcast_indices(As...)
        Mt = imagemeta(As...)
        M = Mt[1]
        if length(Mt) > 1
            for i = 2:length(Mt)
                Mt[i] == M || return ambigop(Symbol(f))
            end
        end
        broadcast!(f, shareproperties(M, similar(data(M), T, shape)), As...)
    end
    # Select all the ImageMeta arrays
    @inline imagemeta(As...) = _imagemeta((), As...)
    _imagemeta(out) = out
    @inline _imagemeta(out, A::ImageMeta, As...) = _imagemeta((out..., A), As...)
    @inline _imagemeta(out, A, As...) = _imagemeta(out, As...)
end

if dotops
    import Base: .+, .-, .*, ./, .^, .<, .>, .==
    batch1 = (:+, :.+, :-, :.-, :*, :.*)
    batch2 = (:+, :.+, :-, :.-)
else
    batch1 = (:+, :-, :*)
    batch2 = (:+, :-)
end

(-)(img::ImageMeta) = shareproperties(img, -data(img))

for op in batch1
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

for op in batch2
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

if dotops
    # evaling in a string avoids a parser depwarn
    include_string("""
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
    """)
end

ambigop(s::Symbol) = error("$s with two ImageMeta arrays: dictionary choice is ambiguous")

using ColorVectorSpace: AbstractGray, TransparentGray, TransparentRGB

# Specialize the low-level broadcasting machinery for ImageMeta.
# This make it work for all operators, rather than specializing
# each operator individually.
Base.BroadcastStyle(::Type{<:ImageMeta}) = Broadcast.ArrayStyle{ImageMeta}()
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{ImageMeta}}, ::Type{ElType}) where ElType
    Mt = imagemeta(bc.args...)
    if length(Mt) > 0
        M = Mt[1]
        if length(Mt) > 1
            for i = 2:length(Mt)
                Mt[i] == M || return ambigop(:Broadcast)
            end
        end
        # Use the properties field of img to create the output
        ret = ImageMeta(similar(Array{ElType}, axes(bc)), M.properties)
    else
        ret = ImageMeta(similar(Array{ElType}, axes(bc)))
    end
    ret
end
# Select all the ImageMeta arrays
@inline imagemeta(As...) = _imagemeta((), As...)
_imagemeta(out) = out
@inline _imagemeta(out, A::ImageMeta, As...) = _imagemeta((out..., A), As...)
@inline _imagemeta(out, A, As...) = _imagemeta(out, As...)

batch1 = (:+, :-, :*)
batch2 = (:+, :-)

(-)(img::ImageMeta) = shareproperties(img, -data(img))

import Base.Broadcast: broadcasted, materialize

for op in batch1
    @eval begin
        broadcasted(::typeof($op),img::ImageMeta{Bool}, n::Bool) = shareproperties(img, materialize(broadcasted(($op),data(img), n)))
        broadcasted(::typeof($op),n::Bool, img::ImageMeta{Bool}) = shareproperties(img, materialize(broadcasted(($op),n, data(img))))
        broadcasted(::typeof($op),img::ImageMeta, n::Number) = shareproperties(img, materialize(broadcasted(($op),data(img), n)))
        broadcasted(::typeof($op),n::Number, img::ImageMeta) = shareproperties(img, materialize(broadcasted(($op),n, data(img))))
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
            ($op)(img::ImageMeta{CV}, n::$CV) where {CV<:$CV} =
                shareproperties(img, ($op)(data(img), n))
            ($op)(n::$CV, img::ImageMeta{CV}) where {CV<:$CV} =
                shareproperties(img, ($op)(n, data(img)))
        end
    end
end

(/)(img::ImageMeta, n::Number) = shareproperties(img, data(img)/n)

ambigop(s::Symbol) = error("$s with two ImageMeta arrays: dictionary choice is ambiguous")

# have to import this or @deprecate doesn't work
import Base: getindex, setindex!, delete!, haskey, get, copy!

@deprecate(ImageMeta(data::AbstractArray{T,N}, props::AbstractDict{String,Any}) where {T,N},
           ImageMeta(data, to_dict(props)))

@deprecate(ImageMeta(data::AbstractArray{T,N}, props::Dict{String,Any}) where {T,N},
           ImageMeta(data, to_dict(props)))

@deprecate(ImageMeta(data::AbstractArray, props::Dict{<:AbstractString}),
           ImageMeta(data, to_dict(props)))

@deprecate(getindex(img::ImageMeta, propname::AbstractString),
           getindex(img, Symbol(propname)))

@deprecate(setindex!(img::ImageMeta, X, propname::AbstractString),
           setindex!(img, X, Symbol(propname)))

@deprecate(copy!(imgdest::ImageMeta, imgsrc::ImageMeta, prop1::AbstractString, props::AbstractString...),
           copy!(imgdest, imgsrc, Symbol(prop1), Symbol.(props)...))

@deprecate(delete!(img::ImageMeta, propname::AbstractString),
           delete!(img, Symbol(propname)))

@deprecate(haskey(img::ImageMeta, k::AbstractString),
           haskey(img, Symbol(k)))

@deprecate(get(img::ImageMeta, k::AbstractString, default),
           get(img, k, Symbol(default)))

function to_dict(dold::AbstractDict{String,Any})
    dnew = Dict{Symbol,Any}()
    for (k, v) in dold
        dnew[Symbol(k)] = v
    end
    return dnew
end


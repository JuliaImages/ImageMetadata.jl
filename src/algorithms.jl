Base.minimum(img::ImageMeta) = minimum(img.data)
Base.maximum(img::ImageMeta) = maximum(img.data)
Base.squeeze(img, dims) = shareproperties(img, squeeze(img.data))

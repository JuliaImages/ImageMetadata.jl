var documenterSearchIndex = {"docs": [

{
    "location": "intro.html#",
    "page": "ImageMetadata.jl",
    "title": "ImageMetadata.jl",
    "category": "page",
    "text": ""
},

{
    "location": "intro.html#ImageMetadata.jl-1",
    "page": "ImageMetadata.jl",
    "title": "ImageMetadata.jl",
    "category": "section",
    "text": "ImageMetadata allows you to add metadata to images: for example, the date and time at which it was collected, identifiers for the location or subject, etc. This metadata is stored as a dictionary, and the ImageMeta type combines properties of arrays and Dict."
},

{
    "location": "intro.html#Introduction-1",
    "page": "ImageMetadata.jl",
    "title": "Introduction",
    "category": "section",
    "text": "You typically create an ImageMeta using keyword arguments:julia> using Colors, ImageMetadata\n\njulia> img = ImageMeta(fill(RGB(1,0,0), 3, 2), date=Date(2016, 7, 31), time=\"high noon\")\nRGB ImageMeta with:\n  data: 3×2 Array{ColorTypes.RGB{FixedPointNumbers.Normed{UInt8,8}},2}\n  properties:\n    time: high noon\n    date: 2016-07-31DocTestSetup = quote\n    using Colors, ImageMetadata\n    img = ImageMeta(fill(RGB(1,0,0), 3, 2), date=Date(2016, 7, 31), time=\"high noon\")\nendYou can then index elements of img like this:julia> img[1,2]\nRGB{N0f8}(1.0,0.0,0.0)and access and set properties like this:julia> img[\"time\"]\n\"high noon\"\n\njulia> img[\"time\"] = \"evening\"\n\"evening\"\n\njulia> img\nRGB ImageMeta with:\n  data: 3×2 Array{ColorTypes.RGB{FixedPointNumbers.Normed{UInt8,8}},2}\n  properties:\n    time: evening\n    date: 2016-07-31You can extract the data matrix with data(img):julia> data(img)\n3×2 Array{ColorTypes.RGB{FixedPointNumbers.Normed{UInt8,8}},2}:\n RGB{N0f8}(1.0,0.0,0.0)  RGB{N0f8}(1.0,0.0,0.0)\n RGB{N0f8}(1.0,0.0,0.0)  RGB{N0f8}(1.0,0.0,0.0)\n RGB{N0f8}(1.0,0.0,0.0)  RGB{N0f8}(1.0,0.0,0.0)and the properties dictionary with properties:julia> properties(img)\nDict{String,Any} with 2 entries:\n  \"time\" => \"high noon\"\n  \"date\" => 2016-07-31Properties are not accessed or modified by most of Images\' algorithms–-the traits that most affect processing are encoded through Julia\'s type system.  However, functions that receive an ImageMeta should return an ImageMeta when appropriate. Naturally, in your own code it\'s fine to use properties to your advantage for custom tasks."
},

{
    "location": "intro.html#Vector-indexing-(region-of-interest-selection)-1",
    "page": "ImageMetadata.jl",
    "title": "Vector indexing (region-of-interest selection)",
    "category": "section",
    "text": "When indexing over an extended area, img[i,j,...] returns an ImageMeta:julia> c = img[1:2, 1:2]\nRGB ImageMeta with:\n  data: 2×2 Array{ColorTypes.RGB{FixedPointNumbers.Normed{UInt8,8}},2}\n  properties:\n    time: high noon\n    date: 2016-07-31This copies both the data (just the relevant portions) and the properties dictionary. In contrast,julia> v = view(img, 1:2, 1:2)\nRGB ImageMeta with:\n  data: 2×2 SubArray{ColorTypes.RGB{FixedPointNumbers.Normed{UInt8,8}},2,Array{ColorTypes.RGB{FixedPointNumbers.Normed{UInt8,8}},2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}\n  properties:\n    time: high noon\n    date: 2016-07-31shares both the data and the properties with the original image img. Modifying values or properties in c has no impact on img, but modifying values or properties in v does."
},

{
    "location": "intro.html#copyproperties/shareproperties-1",
    "page": "ImageMetadata.jl",
    "title": "copyproperties/shareproperties",
    "category": "section",
    "text": "Two convenient ways to construct a new image with the \"same\" properties are copyproperties (makes a copy of the properties dictionary) and shareproperties (shares the properties dictionary).Incidentally, similar makes a copy of the properties dictionary."
},

{
    "location": "intro.html#spatialproperties-1",
    "page": "ImageMetadata.jl",
    "title": "spatialproperties",
    "category": "section",
    "text": "Occasionally you may have a property that is linked to the spatial axes of the image. In such cases, one source for potential confusion is permutedims, which swaps the order of the dimensions in the array: if the order is not also swapped in the appropriate properties, chaos could result.You can declare that certain properties are coupled to spatial axes using \"spatialproperties\":julia> using ImageMetadata\n\njulia> A = reshape(1:15, 3, 5)\n3×5 Base.ReshapedArray{Int64,2,UnitRange{Int64},Tuple{}}:\n 1  4  7  10  13\n 2  5  8  11  14\n 3  6  9  12  15\n\njulia> img = ImageMeta(A, spatialproperties=Set([\"maxsum\"]), maxsum=[maximum(sum(A,1)), maximum(sum(A,2))])\nInt64 ImageMeta with:\n  data: 3×5 Base.ReshapedArray{Int64,2,UnitRange{Int64},Tuple{}}\n  properties:\n    maxsum: [42,45]\n    spatialproperties: Set(String[\"maxsum\"])\n\njulia> imgp = permutedims(img, (2,1))\nInt64 ImageMeta with:\n  data: 5×3 Array{Int64,2}\n  properties:\n    maxsum: [45,42]\n    spatialproperties: Set(String[\"maxsum\"])\n\njulia> maximum(sum(imgp,1))\n45It\'s not possible to anticipate all the possible transformations that might be necessary, but at least simple swaps are handled automatically."
},

{
    "location": "reference.html#",
    "page": "Reference",
    "title": "Reference",
    "category": "page",
    "text": ""
},

{
    "location": "reference.html#ImageMetadata.ImageMeta",
    "page": "Reference",
    "title": "ImageMetadata.ImageMeta",
    "category": "type",
    "text": "ImageMeta is an AbstractArray that can have metadata, stored in a dictionary.\n\nConstruct an image with ImageMeta(A, props) (for a properties dictionary props), or with ImageMeta(A, prop1=val1, prop2=val2, ...).\n\n\n\n\n\n"
},

{
    "location": "reference.html#ImageAxes.data",
    "page": "Reference",
    "title": "ImageAxes.data",
    "category": "function",
    "text": "data(img::ImageMeta) -> array\n\nExtract the data from img, omitting the properties dictionary. array shares storage with img, so changes to one affect the other.\n\nSee also: properties.\n\n\n\n\n\n"
},

{
    "location": "reference.html#ImageMetadata.properties",
    "page": "Reference",
    "title": "ImageMetadata.properties",
    "category": "function",
    "text": "properties(imgmeta) -> props\n\nExtract the properties dictionary props for imgmeta. props shares storage with img, so changes to one affect the other.\n\nSee also: data.\n\n\n\n\n\n"
},

{
    "location": "reference.html#ImageMetadata.copyproperties",
    "page": "Reference",
    "title": "ImageMetadata.copyproperties",
    "category": "function",
    "text": "copyproperties(img::ImageMeta, data) -> imgnew\n\nCreate a new \"image,\" copying the properties dictionary of img but using the data of the AbstractArray data. Note that changing the properties of imgnew does not affect the properties of img.\n\nSee also: shareproperties.\n\n\n\n\n\n"
},

{
    "location": "reference.html#ImageMetadata.shareproperties",
    "page": "Reference",
    "title": "ImageMetadata.shareproperties",
    "category": "function",
    "text": "shareproperties(img::ImageMeta, data) -> imgnew\n\nCreate a new \"image,\" reusing the properties dictionary of img but using the data of the AbstractArray data. The two images have synchronized properties; modifying one also affects the other.\n\nSee also: copyproperties.\n\n\n\n\n\n"
},

{
    "location": "reference.html#ImageMetadata.spatialproperties",
    "page": "Reference",
    "title": "ImageMetadata.spatialproperties",
    "category": "function",
    "text": "spatialproperties(img)\n\nReturn a vector of strings, containing the names of properties that have been declared \"spatial\" and hence should be permuted when calling permutedims.  Declare such properties like this:\n\nimg[\"spatialproperties\"] = [\"spacedirections\"]\n\n\n\n\n\n"
},

{
    "location": "reference.html#ImageCore.spacedirections",
    "page": "Reference",
    "title": "ImageCore.spacedirections",
    "category": "function",
    "text": "spacedirections(img)\n\nUsing ImageMetadata, you can set this property manually. For example, you could indicate that a photograph was taken with the camera tilted 30-degree relative to vertical using\n\nimg[\"spacedirections\"] = ((0.866025,-0.5),(0.5,0.866025))\n\nIf not specified, it will be computed from pixelspacing(img), placing the spacing along the \"diagonal\".  If desired, you can set this property in terms of physical units, and each axis can have distinct units.\n\n\n\n\n\n"
},

{
    "location": "reference.html#Base.permutedims",
    "page": "Reference",
    "title": "Base.permutedims",
    "category": "function",
    "text": "permutedims(img, perm, [spatialprops])\n\nWhen permuting the dimensions of an ImageMeta, you can optionally specify that certain properties are spatial and they will also be permuted. spatialprops defaults to spatialproperties(img).\n\n\n\n\n\n"
},

{
    "location": "reference.html#Reference-1",
    "page": "Reference",
    "title": "Reference",
    "category": "section",
    "text": "ImageMeta\ndata\nproperties\ncopyproperties\nshareproperties\nspatialproperties\nImageMetadata.spacedirections\npermutedims"
},

]}

using ImageMetadata
using Base.Test

@test isempty(detect_ambiguities(ImageMetadata,ImageAxes,ImageCore,IndirectArrays,Base,Core))

include("core.jl")
include("operations.jl")
info("Beginning of tests with deprecation warnings\n\n")
include("deprecated.jl")

nothing

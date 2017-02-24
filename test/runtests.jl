using ImageMetadata
using Base.Test

if VERSION < v"0.6.0-dev"
    @test isempty(detect_ambiguities(ImageMetadata,ImageAxes,ImageCore,IndirectArrays,Base,Core))
end

include("core.jl")
include("operations.jl")
info("Beginning of tests with deprecation warnings\n\n")
include("deprecated.jl")

nothing

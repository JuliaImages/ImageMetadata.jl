using ImagesMeta
using Base.Test

# write your own tests here
@test isempty(detect_ambiguities(ImagesMeta,ImagesAxes,ImagesCore,IndirectArrays,Base,Core))

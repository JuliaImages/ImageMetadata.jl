using Colors, ColorVectorSpace, SimpleTraits, ImagesAxes, ImagesMeta
using Base.Test

@test isempty(detect_ambiguities(ImagesMeta,ImagesAxes,ImagesCore,IndirectArrays,Base,Core))

@traitimpl TimeAxis{Axis{:time}}

@testset "indexing" begin
    for A in (rand(3,5),
              view(rand(4,6), 1:3, 1:5),
              view(rand(4,5), [1,2,4], :),
              reshape(1:15, 3, 5),
              rand(RGB{U8}, 3, 5),
              rand(Gray{Float32}, 3, 5))
        img = ImageMeta(A; prop1 = 1, prop2 = [1,2,3])
        @test eltype(img) == eltype(A)
        @test ndims(img) == 2
        @test size(img) == (3,5)
        @test data(img) === A
        for j = 1:5, i = 1:3
            @test img[i,j] === A[i,j]
        end
        for k = 1:15
            @test img[k] === A[k]
        end
        k = 0
        for a in img
            @test a == A[k+=1]
        end
        if !isa(A, typeof(reshape(1:15, 3, 5)))
            img[2,3] = zero(eltype(img))
            @test A[2,3] == zero(eltype(A))
            img[4] = one(eltype(img))
            @test A[4] == one(eltype(A))
        end
        @test_throws BoundsError img[0,0]
        @test_throws BoundsError img[4,1]
        @test_throws BoundsError img[1,6]
        @test img["prop1"] == 1
        @test img["prop2"] == [1,2,3]
        img["prop1"] = -1
        @test img["prop1"] == -1
    end
    a = zeros(3)
    sizehint!(a, 10)
    @test_throws BoundsError a[5]
    @inbounds a[5] = 1.234
    @inbounds val = a[5]
    @test val == 1.234
    a = zeros(3,5)
end

@testset "copy/similar" begin
    img = ImageMeta(rand(3,5); prop1 = 1, prop2 = [1,2,3])
    img2 = copy(img)
    @test img2.data == img.data
    img2[2,2] = -1
    @test img2[2,2] < 0
    @test img[2,2] >= 0
    img2["prop2"][2] = -2
    @test img2["prop2"] == [1,-2,3]
    @test img["prop2"] == [1,2,3]
    img2 = similar(img)
    @test img2["prop1"] == 1
    @test img2["prop2"] == [1,2,3]
    @test img2.data != img.data
    img2["prop3"] = 7
    @test !haskey(img, "prop3")
end

@testset "copy/shareproperties" begin
    img = ImageMeta(rand(3,5); prop1 = 1, prop2 = [1,2,3])
    @test !isempty(properties(img))
    img2 = copyproperties(img, reshape(1:15, 5, 3))
    @test size(img2) == (5,3)
    img2["prop1"] = -1
    @test img["prop1"] == 1
    img2 = shareproperties(img, reshape(1:15, 5, 3))
    @test size(img2) == (5,3)
    img2["prop1"] = -1
    @test img["prop1"] == -1
    imgb = ImageMeta(rand(RGB{U8}, 2, 2), propa = "hello", propb = [1,2])
    copy!(img, imgb, "propa", "propb")
    @test img["propa"] == "hello"
    @test img["propb"] == [1,2]
    img["propb"][2] = 10
    @test img["propb"] == [1,10]
    @test imgb["propb"] == [1,2]
    delete!(img, "propb")
    @test  haskey(img, "propa")
    @test !haskey(img, "propb")
end

@testset "show" begin
    for supp in (Set(["prop3"]), "prop3")
        img = ImageMeta(rand(3,5); prop1 = 1, prop2 = [1,2,3], suppress = supp, prop3 = "hide")
        io = IOBuffer()
        show(io, img)
        str = takebuf_string(io)
        @test contains(str, "ImageMeta with")
        @test contains(str, "prop1")
        @test contains(str, "prop2")
        @test contains(str, "[1,2,3]")
        @test contains(str, "prop3")
        @test !contains(str, "hide")
        @test contains(str, "<suppressed>")
        @test !contains(str, "suppress:")
    end
end

@testset "traits" begin
    img = ImageMeta(AxisArray(rand(3,5,8),
                              Axis{:x}(1:3),
                              Axis{:y}(1:5),
                              Axis{:time}(0.1:0.1:0.8)))
    @test @inferred(timedim(img)) == 3
    @test @inferred(pixelspacing(img)) === (1,1,0.1)
    @test spacedirections(img) === ((1,0),(0,1))
    img["spacedirections"] = "broken"
    @test spacedirections(img) == "broken"
    @test sdims(img) == 2
    @test @inferred(coords_spatial(img)) == (1,2)
    @test spatialorder(img) == (:x, :y)
    @test nimages(img) == 8
    @test @inferred(size_spatial(img)) == (3,5)
    @test @inferred(indices_spatial(img)) == (Base.OneTo(3), Base.OneTo(5))
    assert_timedim_last(img)
end

include("deprecated.jl")

nothing

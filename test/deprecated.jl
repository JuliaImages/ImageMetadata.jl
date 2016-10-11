using FixedPointNumbers, Colors, ColorVectorSpace, SimpleTraits, ImageAxes, ImageMetadata
using Base.Test

msg_contains(pass, msg) = contains(pass.value.msg, msg) || error(pass.value.msg, " does not contain \"", msg, "\"")

@testset "Deprecated" begin
    a = rand(3,3)
    @inferred(Image(a))
    B = rand(convert(UInt16, 1):convert(UInt16, 20), 3, 5)
    cmap = reinterpret(RGB, repmat(reinterpret(N0f8, round(UInt8, linspace(12, 255, 20)))', 3, 1))
    img = ImageCmap(copy(B), cmap, prop="cmap")
    imgd = copy(img)
    @test size(img) == size(imgd) == (3,5)
    @test ndims(img) == ndims(imgd) == 2

    @testset "indexing" begin
        c = img[1,1]
        @test red(c) == green(c) == blue(c)
        @test_throws ErrorException (img[1,1] = RGB{N0f8}(0.2,0.4,0.6))
        imgd[1,1] = RGB{N0f8}(0.2,0.4,0.6)
        @test imgd[1,1] == RGB{N0f8}(0.2,0.4,0.6)
    end

    @testset "views" begin
        sl = subim(img, 1:2, 1:2)
        @test haskey(sl, "prop")
        sl = sliceim(img, 1:2, 1:2)
        @test haskey(sl, "prop")
        imga = ImageMeta(AxisArray(rand(3,5,8),
                                   Axis{:x}(1:3),
                                   Axis{:y}(1:5),
                                   Axis{:time}(0.1:0.1:0.8)),
                         dummy=true)
        @test view(imga, "y", 2:3) == imga.data[Axis{:y}(2:3)]
        v = viewim(imga, "y", 2)
        @test v["dummy"] == true
        v = subim(imga, "y", 2)
        @test v["dummy"] == true
        v = sliceim(imga, "y", 2)
        @test v["dummy"] == true
    end

    @testset "traits" begin
        @test  isdirect(imgd)
        @test !isdirect(img)
        @test colorspace(B) == "Gray"
        @test colorspace(img) == "RGB"
        @test colorspace(imgd) == "RGB"
        @test ncolorelem(img) == ncolorelem(imgd) == 1
        @test nimages(img) == nimages(imgd) == 1
    end

    @testset "deprecated properties" begin
        @test isempty(properties(rand(3,3)))
        for Im in (Image, ImageMeta)
            result = @test_throws ErrorException Im(rand(3,5,5), colorspace="RGB")
            msg_contains(result, "color is encoded")
            Im(rand(3,5,5), clrspace="RGB")  # ensure it's specific for that name
            result = @test_throws ErrorException Im(rand(3,5,5), colordim=1)
            msg_contains(result, "color is encoded")
            result = @test_throws ErrorException Im(rand(3,5), limits=(0.25,0.75))
            msg_contains(result, "limits are always")
            result = @test_throws ErrorException Im([1:3, 4:7], limits=(0.25,0.75)) # case where zero(eltype) fails
            msg_contains(result, "limits are always")
            result = @test_throws ErrorException Im(rand(3,5), pixelspacing=[2,1])
            msg_contains(result, "please switch to ImageAxes")
            result = @test_throws ErrorException Im(rand(3,5,12), timedim=3)
            msg_contains(result, "please switch to ImageAxes")
            result = @test_throws ErrorException Im(rand(3,5), spatialorder=["boo", "rah"])
            msg_contains(result, "data, :boo, :rah")
        end
    end
end

nothing

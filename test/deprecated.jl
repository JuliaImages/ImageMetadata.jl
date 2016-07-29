using Colors, ColorVectorSpace, SimpleTraits, ImagesAxes, ImagesMeta
using Base.Test

testing_units = Int == Int64
if testing_units
    using SIUnits, SIUnits.ShortUnits
end

@testset "Deprecated" begin
    a = rand(3,3)
    @inferred(Image(a))
    B = rand(convert(UInt16, 1):convert(UInt16, 20), 3, 5)
    cmap = reinterpret(RGB, repmat(reinterpret(U8, round(UInt8, linspace(12, 255, 20)))', 3, 1))
    img = ImageCmap(copy(B), cmap, Dict{String, Any}("pixelspacing"=>[2.0, 3.0],
                                                     "spatialorder"=>ImagesCore.yx))
    imgd = copy(img)

    @testset "indexing" begin
        c = img[1,1]
        @test red(c) == green(c) == blue(c)
        @test_throws ErrorException (img[1,1] = RGB{U8}(0.2,0.4,0.6))
        imgd[1,1] = RGB{U8}(0.2,0.4,0.6)
        @test imgd[1,1] == RGB{U8}(0.2,0.4,0.6)
    end


    @test size(img) == size(imgd) == (3,5)
    @test ndims(img) == ndims(imgd) == 2

    @testset "traits" begin
        @test  isdirect(imgd)
        @test !isdirect(img)
        @test colorspace(B) == "Gray"
        @test colorspace(img) == "RGB"
        @test colorspace(imgd) == "RGB"
        @test ncolorelem(img) == ncolorelem(imgd) == 3
        @test nimages(img) == nimages(imgd) == 1
    end

    @testset "deprecated properties" begin
        f, io = mktemp()
        OLDSTDERR = STDERR
        redirect_stderr(io)
        @test pixelspacing(img) == (1,1)
        redirect_stderr(OLDSTDERR)
        close(io)
        str = readstring(f)
        @show str
        @test contains(str, "pixelspacing property")
        if testing_units
            imgd["pixelspacing"] = [2.0mm, 3.0mm]
            @test imgd["pixelspacing"] == [2.0mm,3.0mm]
        end
        A = AxisArray(rand(3,5), Axis{:x}(1mm:2mm:5mm), Axis{:y}(1mm:1mm:5mm))
        imgax = Image(A, pixelspacing=[0,0])
        @test pixelspacing(imgax) == (2mm,1mm)
    end

end

nothing

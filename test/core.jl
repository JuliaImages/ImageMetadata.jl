using FixedPointNumbers, Colors, ColorVectorSpace, SimpleTraits, ImageAxes, ImageMetadata
using Base.Test

@testset "indexing" begin
    # 1d images
    for A = (rand(3),
             view(rand(6), 1:2:5),
             view(rand(6), [1,2,4]),
             rand(RGB{N0f8}, 3),
             rand(Gray{N0f8}, 3))
        img = ImageMeta(A, prop1=1, prop2=[1,2,3])
        @test eltype(img) == eltype(A)
        @test Base.linearindexing(img) == Base.linearindexing(A)
        @test ndims(img) == 1
        @test size(img) == (3,)
        @test data(img) === A
        for i = 1:3
            @test img[i] === A[i]
        end
        for I in eachindex(img)
            @test img[I] === A[I]
        end
        k = 0
        for a in img
            @test a == A[k+=1]
        end
        img[2] = zero(eltype(img))
        @test A[2] == zero(eltype(A))
        img[3] = one(eltype(img))
        @test A[3] == one(eltype(A))
        @test_throws BoundsError img[0]
        @test_throws BoundsError img[4]
        @test img["prop1"] == 1
        @test img["prop2"] == [1,2,3]
        img["prop1"] = -1
        @test img["prop1"] == -1
    end

    for A in (rand(3,5),
              view(rand(4,6), 1:3, 1:5),
              view(rand(4,5), [1,2,4], :),
              reshape(1:15, 3, 5),
              rand(RGB{N0f8}, 3, 5),
              rand(Gray{Float32}, 3, 5))
        img = ImageMeta(A, prop1=1, prop2=[1,2,3]) # TODO: add @inferred (see julia #17719)
        @test eltype(img) == eltype(A)
        @test Base.linearindexing(img) == Base.linearindexing(A)
        @test ndims(img) == 2
        @test size(img) == (3,5)
        @test data(img) === A
        for j = 1:5, i = 1:3
            @test img[i,j] === A[i,j]
        end
        for k = 1:15
            @test img[k] === A[k]
        end
        for I in eachindex(img)
            @test img[I] === A[I]
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
    # Test bounds-checking removal by @inbounds
    if Base.JLOptions().check_bounds != 1 && Base.JLOptions().can_inline == 1
        a = zeros(3)
        sizehint!(a, 10)  # make sure we don't cause a segfault
        @test_throws BoundsError a[5]
        @inbounds a[5] = 1.234
        @inbounds val = a[5]
        @test val == 1.234
        a = zeros(3,5)
    end

    # Issue #10
    A = AxisArray(rand(3,5), :y, :x)
    B = ImageMeta(A, info="blah")
    Broi = B[2:3, 2:3]
    @test isa(Broi, AxisArray)
    @test axisnames(Broi) == (:y, :x)
    A1, B1 = A[2:7], B[2:7]
    @test typeof(A1) == typeof(B1) && A1 == B1
    Broi = view(B, 2:3, 2:3)
    @test isa(Broi, AxisArray)
    @test axisnames(Broi) == (:y, :x)
end

@testset "convert" begin
    A = rand(3,5)
    M = convert(ImageMeta, A)
    @test isa(M, ImageMeta)
    @test convert(ImageMeta, M) === M
    @test convert(ImageMeta{Float64}, M) === M
    @test eltype(convert(ImageMeta{Gray{Float64}}, M)) == Gray{Float64}
    @test eltype(convert(ImageMeta{Gray}, M)) == Gray{Float64}
    @test convert(ImageMeta{Gray}, M) == M
    @test convert(ImageMeta{Gray}, A) == A
end

@testset "reinterpret" begin
    # It's possible that reinterpret shouldn't be defined for ImageMeta, but...
    A = rand(Float32, 4, 5)
    M = ImageMeta(A, meta=true)
    Mr = reinterpret(Gray, M)
    @test eltype(Mr) == Gray{Float32}
    @test Mr["meta"] = true
    # Ensure that it never gets defined for the un-reinterpretable
    M = ImageMeta(view(A, 1:2:3, 1:4), meta=true)
    @test_throws ErrorException reinterpret(Gray, M)
    @test_throws ErrorException reinterpret(Gray{Float32}, M)
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
    img2 = similar(img, RGB{Float16}, (Base.OneTo(5),))
    @test img2["prop1"] == 1
    @test img2["prop2"] == [1,2,3]
    @test eltype(img2) == RGB{Float16}
    @test size(img2) == (5,)
    @test img2.data != img.data
    img2["prop3"] = 7
    @test !haskey(img, "prop3")

    A = AxisArray(rand(3,5), :y, :x)
    B = ImageMeta(A, info="blah")
    C = similar(B)
    @test isa(C, ImageMeta)
    @test isa(C.data, AxisArray)
    @test eltype(C) == Float64
    C = similar(B, RGB{Float16})
    @test isa(C, ImageMeta)
    @test isa(C.data, AxisArray)
    @test eltype(C) == RGB{Float16}
end

@testset "copy/shareproperties/viewim" begin
    img = ImageMeta(rand(3,5); prop1 = 1, prop2 = [1,2,3])
    @test !isempty(properties(img))
    v = viewim(img, 1:2, 1:2)
    c = getindexim(img, 1:2, 1:2)
    @test v["prop1"] == 1
    @test c["prop1"] == 1
    img2 = copyproperties(img, reshape(1:15, 5, 3))
    @test size(img2) == (5,3)
    img2["prop1"] = -1
    @test img["prop1"] == 1
    img2 = shareproperties(img, reshape(1:15, 5, 3))
    @test size(img2) == (5,3)
    img2["prop1"] = -1
    @test img["prop1"] == -1
    @test v["prop1"] == -1
    @test c["prop1"] == 1
    imgb = ImageMeta(rand(RGB{N0f8}, 2, 2), propa = "hello", propb = [1,2])
    copy!(img, imgb, "propa", "propb")
    @test img["propa"] == "hello"
    @test img["propb"] == [1,2]
    img["propb"][2] = 10
    @test img["propb"] == [1,10]
    @test imgb["propb"] == [1,2]
    delete!(img, "propb")
    @test  haskey(img, "propa")
    @test !haskey(img, "propb")

    img = ImageMeta(AxisArray(rand(3,5,8),
                              Axis{:x}(1:3),
                              Axis{:y}(1:5),
                              Axis{:time}(0.1:0.1:0.8));
                    prop1 = 1, prop2 = [1,2,3])
    v = viewim(img, Axis{:time}(0.25..0.5))
    c = getindexim(img, Axis{:time}(0.25..0.5))
    @test v["prop1"] == 1
    @test c["prop1"] == 1
end

@testset "views" begin
    for A in (rand(Gray{N0f8}, 4, 5), rand(RGB{Float32}, 4, 5))
        t = now()
        M = ImageMeta(A, date=t)
        vM = channelview(M)
        @test isa(vM, ImageMeta)
        @test vM["date"] == t
        @test data(vM) == channelview(A)
        @test isa(rawview(vM), ImageMeta)
        @test rawview(vM)[1,2] === rawview(channelview(A))[1,2]
    end
    for (A,C) in ((rand(UInt8, 4, 5), Gray), (rand(UInt8, 3, 4, 5), RGB))
        M = ImageMeta(A)
        vM = normedview(M)
        @test isa(vM, ImageMeta)
        @test data(vM) == normedview(A)
        cvM = colorview(C, vM)
        @test isa(cvM, ImageMeta)
        @test cvM[1,2] === colorview(C, normedview(A))[1,2]
    end
    A = AxisArray(rand(RGB{N0f8}, 3, 5), :y, :x)
    t = now()
    M = ImageMeta(A, date=t)
    vM = channelview(M)
    @test isa(vM, ImageMeta)
    @test colordim(vM) == 1
    pvM = permutedims(vM, (2,3,1))
    @test colordim(pvM) == 3
    @test pvM["date"] == t
    sleep(0.1)
    pvM["date"] = now()
    @test M["date"] == t && pvM["date"] != t
    vM = permuteddimsview(M, (2,1))
    @test vM["date"] == t
    @test axisnames(vM) == (:x, :y)
end

@testset "meta-axes" begin
    A = AxisArray(rand(3,5), :y, :x)
    M = ImageMeta(A)
    @test axes(M) == (Axis{:y}(1:3), Axis{:x}(1:5))
    @test axisdim(M, Axis{:y}) == 1
    @test axisdim(M, Axis{:x}) == 2
    @test_throws ErrorException axisdim(M, Axis{:z})
    @test axisnames(M) == (:y, :x)
    @test axisvalues(M) == (1:3, 1:5)
    @test M[Axis{:y}(2:3)] == A[Axis{:y}(2:3)]
    @test view(M, Axis{:y}(2:3)) == A[Axis{:y}(2:3)]
    M[Axis{:y}(2:3), Axis{:x}(1)] = -5
    @test all(A[2:3,1] .== -5)
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
        io = IOBuffer()
        show(io, MIME("text/plain"), img)
        str2 = takebuf_string(io)
        @test str == str2
    end
end

@testset "traits" begin
    img = ImageMeta(AxisArray(rand(3,5,8),
                              Axis{:x}(1:3),
                              Axis{:y}(1:5),
                              Axis{:time}(0.1:0.1:0.8)))
    @test @inferred(timedim(img)) == 3
    @test @inferred(pixelspacing(img)) === (1,1)
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

@testset "spatialprops" begin
    img = ImageMeta(rand(3,5),
                    spatialproperties=Set(["vector","matrix","tuple"]),
                    vector=[1,2],
                    matrix=[1 3; 2 4],
                    tuple=(1,2))
    for imgp in (img', permutedims(img, (2,1)), permuteddimsview(img, (2,1)))
        @test imgp.data == img.data'
        @test imgp["vector"] == [2,1]
        @test imgp["matrix"] == [4 2; 3 1]
        @test imgp["tuple"]  == (2,1)
        @test img["vector"] == [1,2]
        @test img["matrix"] == [1 3; 2 4]
        @test img["tuple"]  == (1,2)
    end
end

@testset "dimchange" begin
    M = ImageMeta([1,2,3,4])
    Mp = M'
    @test ndims(Mp) == 2 && isa(Mp, ImageMeta)

    M = ImageMeta([1,2,3,4],
                  spatialproperties=["vector"],
                  vector=[1])
    @test_throws ErrorException M'
end

nothing

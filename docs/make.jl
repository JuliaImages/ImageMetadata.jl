using Documenter, ImageMetadata

makedocs(modules  = [ImageMetadata],
         format   = Documenter.Formats.HTML,
         sitename = "ImageMetadata",
         pages    = ["intro.md", "reference.md"])

deploydocs(repo   = "github.com/JuliaImages/ImageMetadata.jl.git",
           julia  = "0.5",
           target = "build",
           deps   = nothing,
           make   = nothing)
#           deps   = Deps.pip("mkdocs", "python-markdown-math"),

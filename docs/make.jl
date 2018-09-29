using Documenter, ImageMetadata

makedocs(modules  = [ImageMetadata],
         format   = :html,
         sitename = "ImageMetadata",
         pages    = ["intro.md", "reference.md"])

deploydocs(repo   = "github.com/JuliaImages/ImageMetadata.jl.git",
           julia  = "1.0",
           target = "build",
           deps   = nothing,
           make   = nothing)
#           deps   = Deps.pip("mkdocs", "python-markdown-math"),

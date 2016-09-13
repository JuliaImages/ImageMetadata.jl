using Documenter, ImageMetadata

makedocs(modules = [ImageMetadata])

deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
           repo   = "github.com/JuliaImages/ImageMetadata.jl.git")

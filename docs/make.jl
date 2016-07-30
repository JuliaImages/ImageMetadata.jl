using Documenter, ImagesMeta

makedocs(modules = [ImagesMeta])

deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
           repo   = "github.com/JuliaImages/ImagesMeta.jl.git")

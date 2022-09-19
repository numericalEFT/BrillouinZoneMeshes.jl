using BZMeshes
using Documenter

DocMeta.setdocmeta!(BZMeshes, :DocTestSetup, :(using BZMeshes); recursive=true)

makedocs(;
    modules=[BZMeshes],
    authors="Xiansheng Cai, Tao Wang and Kun Chen",
    repo="https://github.com/numericaleft/BZMeshes.jl/blob/{commit}{path}#{line}",
    sitename="BZMeshes.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://numericaleft.github.io/BZMeshes.jl",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(;
    repo="github.com/numericalEFT/BZMeshes.jl",
    devbranch="master"
)

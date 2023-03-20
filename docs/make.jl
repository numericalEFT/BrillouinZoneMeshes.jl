using BrillouinZoneMeshes
using Documenter

DocMeta.setdocmeta!(BrillouinZoneMeshes, :DocTestSetup, :(using BrillouinZoneMeshes); recursive=true)

makedocs(;
    modules=[BrillouinZoneMeshes],
    authors="Xiansheng Cai, Tao Wang and Kun Chen",
    repo="https://github.com/numericaleft/BrillouinZoneMeshes.jl/blob/{commit}{path}#{line}",
    sitename="BrillouinZoneMeshes.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://numericaleft.github.io/BrillouinZoneMeshes.jl",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => "api_reference.md",
        "SubModules" => [
            "AbstractMeshes" => "lib/AbstractMeshes.md",
            "Cells" => "lib/Cells.md",
            "BaseMesh" => "lib/BaseMesh.md",
            "MeshMaps" => "lib/MeshMaps.md",
            "BZMeshes" => "lib/BZMeshes.md",
        ]
    ]
)

deploydocs(;
    repo="github.com/numericalEFT/BrillouinZoneMeshes.jl",
    devbranch="master"
)

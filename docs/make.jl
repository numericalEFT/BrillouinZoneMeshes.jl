using SpaceGrid
using Documenter

DocMeta.setdocmeta!(SpaceGrid, :DocTestSetup, :(using SpaceGrid); recursive=true)

makedocs(;
    modules=[SpaceGrid],
    authors="Tao Wang, Xiansheng Cai",
    repo="https://github.com/fsxbhyy/SpaceGrid.jl/blob/{commit}{path}#{line}",
    sitename="SpaceGrid.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fsxbhyy.github.io/SpaceGrid.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fsxbhyy/SpaceGrid.jl",
    devbranch="main",
)

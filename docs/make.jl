using ImplicitIntegration
using Documenter

DocMeta.setdocmeta!(ImplicitIntegration, :DocTestSetup, :(using ImplicitIntegration); recursive=true)

makedocs(;
    modules=[ImplicitIntegration],
    authors="Luiz M. Faria",
    sitename="ImplicitIntegration.jl",
    format=Documenter.HTML(;
        canonical="https://maltezfaria.github.io/ImplicitIntegration.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/maltezfaria/ImplicitIntegration.jl",
    devbranch="main",
)

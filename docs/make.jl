using ImplicitIntegration
using Documenter
using DocumenterCitations
using GLMakie

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :numeric)

DocMeta.setdocmeta!(
    ImplicitIntegration,
    :DocTestSetup,
    :(using ImplicitIntegration);
    recursive = true,
)

makedocs(;
    modules = [ImplicitIntegration],
    sitename = "ImplicitIntegration.jl",
    format = Documenter.HTML(;
        canonical = "https://maltezfaria.github.io/ImplicitIntegration.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md", "Docstrings" => "docstrings.md"],
    plugins = [bib],
)

deploydocs(; repo = "github.com/maltezfaria/ImplicitIntegration.jl", devbranch = "main")

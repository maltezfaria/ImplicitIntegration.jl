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
    authors = "Luiz M. Faria and Antoine Levitt",
    format = Documenter.HTML(;
        canonical = "https://maltezfaria.github.io/ImplicitIntegration.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Specializations" => "specializations.md",
        "Docstrings" => "docstrings.md",
        "Bibliography" => "bibliography.md",
    ],
    plugins = [bib],
    draft = false,
)

deploydocs(;
    repo = "github.com/maltezfaria/ImplicitIntegration.jl",
    devbranch = "main",
    push_preview = true,
)

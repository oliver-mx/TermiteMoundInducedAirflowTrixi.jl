using TermiteMoundInducedAirflowTrixi
using Documenter

DocMeta.setdocmeta!(
    TermiteMoundInducedAirflowTrixi,
    :DocTestSetup,
    :(using TermiteMoundInducedAirflowTrixi);
    recursive = true,
)

makedocs(;
    modules = [TermiteMoundInducedAirflowTrixi],
    authors = "Oliver Peter Marx <oliver-mx@uni-hamburg.de>",
    sitename = "TermiteMoundInducedAirflowTrixi.jl",
    format = Documenter.HTML(;
        canonical = "https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tuto.md",
        "Mathematics" => "math.md",
        "Authors" => "auth.md",
        "License" => "lice.md",
        "References" => "reference.md",
    ],
)

deploydocs(;
    repo = "github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl",
    devbranch = "main",
)

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
        edit_link = "master",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(;
    repo = "github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl",
    devbranch = "master",
)
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
    authors = "oliver-mx <oliver-mx@uni-hamburg.de>",
    sitename = "TermiteMoundInducedAirflowTrixi.jl",
    format = Documenter.HTML(;
        canonical = "https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl",
        edit_link = "master",
        assets = String[],
    ),
    pages=["Home" => "index.md", "Tutorial" => "tuto.md", "Mathematics" => "math.md"],
)

deploydocs(;
    repo = "github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl",
    devbranch = "master",
)
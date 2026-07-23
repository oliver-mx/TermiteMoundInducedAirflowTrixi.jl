using TermiteMoundInducedAirflowTrixi
using Documenter
using DocumenterInterLinks
using Literate
using Changelog: Changelog

# Provide external links (project root and inventory file)
links = InterLinks("TermiteMoundInducedAirflowTrixi" => ("https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl/"))

DocMeta.setdocmeta!(
    TermiteMoundInducedAirflowTrixi,
    :DocTestSetup,
    :(using TermiteMoundInducedAirflowTrixi);
    recursive = true,
)

# Copy contents form README to the starting page to not need to synchronize it manually
readme_text = read(joinpath(dirname(@__DIR__), "README.md"), String)
readme_text = replace(readme_text,
                      "[LICENSE.md](LICENSE.md)" => "[License](@ref)",
                      "<p" => "```@raw html\n<p",
                      "p>" => "p>\n```",
                      r"\[comment\].*\n" => "")    # remove comments
write(joinpath(@__DIR__, "src", "home.md"), readme_text)

makedocs(;
    modules = [TermiteMoundInducedAirflowTrixi],
    authors = "Oliver Peter Marx <oliver-mx@uni-hamburg.de>",
    sitename = "TermiteMoundInducedAirflowTrixi.jl",
    format = Documenter.HTML(;
        canonical = "https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl",
        edit_link = "master",
        assets = String[],
    ),
    pages = ["Home" => "index.md", 
        "Installation" => "installation.md",
        "Introduction" => "introduction.md",
        "License" => "license.md",
        "Reference" => "reference.md"],
    plugins = [links],)

deploydocs(;
    repo = "github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl",
    devbranch = "main",
)

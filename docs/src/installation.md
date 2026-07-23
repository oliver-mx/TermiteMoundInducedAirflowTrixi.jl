## Installation
If you have not yet installed Julia, please [follow the instructions for your
operating system](https://julialang.org/downloads/platform/). TermiteMoundInducedAirflowTrixi.jl works
with Julia v1.12.6 and newer. We recommend using the latest stable release of Julia.

TermiteMoundInducedAirflowTrixi.jl and its related tools are registered Julia packages. 
Below is a list of relevant packages and sub-packages:

* [Trixi.jl](https://github.com/trixi-framework/Trixi.jl)
* [Trixi2Vtk.jl](https://github.com/trixi-framework/Trixi2Vtk.jl)
* [OrdinaryDiffEqLowStorageRK.jl](https://github.com/SciML/OrdinaryDiffEq.jl)
* [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)
* [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl)
* [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl)
* [Plots.jl](https://github.com/JuliaPlots/Plots.jl)

They can be added and installed by executing the following commands in the Julia REPL:
```julia
julia> using Pkg

julia> Pkg.add(["TermiteMoundInducedAirflowTrixi", "Trixi", "Trixi2Vtk", "OrdinaryDiffEqLowStorageRK",
                "Interpolations", "QuadGK", "FastGaussQuadrature", "Plots"])
```
You can copy and paste all commands to the REPL *including* the leading
`julia>` prompts - they will automatically be stripped away by Julia.

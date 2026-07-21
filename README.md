# TermiteMoundInducedAirflowTrixi

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl/dev/)
[![Build Status](https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/oliver-mx/TermiteMoundInducedAirflowTrixi.jl/graph/badge.svg?token=w06WNAQpNI)](https://codecov.io/gh/oliver-mx/TermiteMoundInducedAirflowTrixi.jl)

**TermiteMoundInducedAirflowTrixi.jl** is a numerical simulation package focused on solving a thermo fluid dynamic model describing airflow inside of a termite mound with the discontinuous Galerkin method written in Julia. The package builds on the numerical simulation framework for conservation laws [Trixi.jl](https://github.com/trixi-framework/Trixi.jl). It provides specialized models for termite mound and passive house applications.

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

## Referencing
You can directly refer to TermiteMoundInducedAirflowTrixi.jl as

```bibtex
@misc{marx26,
  title={{TermiteMoundInducedAirflowTrixi.jl}: {T}ermite mound and passive house airflow simulations with {T}rixi.jl},
  author={Marx, Oliver P and Gasser, Ingenuin and Schmidgall, Annika},
  year={2026},
  howpublished={\url{https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl}},
  doi={https://doi.org/... }
}
```

## Authors
TermiteMoundInducedAirflowTrixi.jl is maintained by [Oliver P. Marx](https://github.com/oliver-mx)

## License
TermiteMoundInducedAirflowTrixi.jl is licensed under the MIT license (see [LICENSE.md](LICENSE.md)).

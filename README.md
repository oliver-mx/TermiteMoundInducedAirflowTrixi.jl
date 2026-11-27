# TermiteMoundInducedAirflowTrixi

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl/dev/)
[![Build Status](https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl/actions/workflows/CI.yml/badge.svg?)](https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/oliver-mx/TermiteMoundInducedAirflowTrixi.jl/graph/badge.svg?token=w06WNAQpNI)](https://codecov.io/gh/oliver-mx/TermiteMoundInducedAirflowTrixi.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22708420.svg)](https://doi.org/10.5281/zenodo.22708420)

**TermiteMoundInducedAirflowTrixi.jl** is a numerical simulation package written in Julia. It solves a thermo-fluid-dynamic model of the airflow inside a termite mound with the discontinuous Galerkin method. The package builds on the numerical simulation framework for conservation laws [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).

## Description

This package provides:

- **1D Termite Mound Model**: This thermo-fluid-dynamic model describes the airflow dynamics inside termite mounds. The one-dimensional mathematical model originates from the Euler equations of gas dynamics in a low Mach number regime:

```math
\begin{aligned}
\frac{\partial \rho}{\partial t} + \frac{\partial (\rho u)}{\partial x}
  &= -\frac{A_x}{A} \rho u, \\
\frac{\partial T_i}{\partial t} &= k_i (T - T_i), \\
\frac{\partial p_0}{\partial t} + \gamma p_0 \frac{\partial u}{\partial x}
  &= -\gamma \frac{A_x}{A} u p_0 - \frac{k_w}{\sqrt{A}} (T - T_u), \\
\frac{\mathrm{d} v}{\mathrm{d} t}
  &= \frac{1}{\int_0^1 \frac{\rho}{A} \, dy} \int_0^1 \Big[ -\rho u u_x
     - \frac{\beta}{A} \rho u \big(\eta + (1-\eta)|u|\big)
     - \frac{h_x}{Fr^2} (\rho - \rho_{h_0}) \Big] \, dy, \\
p_0 &= \rho T, \\
u(t, x) &= \frac{v(t)}{A(x)} + \frac{1}{A(x)} \int_0^x \frac{A(y)}{\gamma p_0}
  \left[ -\frac{\mathrm{d} p_0}{\mathrm{d} t} - \frac{k_w}{\sqrt{A}} (T - T_u) \right] dy.
\end{aligned}
```

- **The equations** `TermiteMoundEquations1D` for [Trixi.jl](https://github.com/trixi-framework/Trixi.jl),
  a high-performance framework for conservation laws with the discontinuous Galerkin (DG)
  method, together with the model parameters (`termite_parameters`) and the model data
  (cross section, height profile and boundary temperature) of termite mounds.
- **Two numerical methods** for the PDE-ODE system with algebraic constraint:
  - a **coupled approach** (`CoupledSourceTerms` + `AlgebraicVelocityCallback`): all
    equations are integrated together. **This is the recommended approach.**
  - a **Godunov-type splitting** (`GodunovSplittingCallback`): a first-order operator
    splitting.
- **A manufactured solution** to verify both methods (see the
  [documentation](https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl/stable/manufactured/)).
- **Example elixirs** in `examples/` for a termite mound, a reference mound and the
  manufactured solution, each for both methods.

To get started, follow the
[tutorial](https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl/stable/tuto/).

## Installation
If you have not yet installed Julia, please [follow the instructions for your
operating system](https://julialang.org/downloads/platform/). TermiteMoundInducedAirflowTrixi.jl works
with Julia v1.10 and newer. We recommend using the latest stable release of Julia.

TermiteMoundInducedAirflowTrixi.jl and its related tools are registered Julia packages.
To run the examples in `examples/`, the following packages are required:

* [Trixi.jl](https://github.com/trixi-framework/Trixi.jl): DG solver framework (mesh, semidiscretization, callbacks)
* [OrdinaryDiffEqLowStorageRK.jl](https://github.com/SciML/OrdinaryDiffEq.jl): low-storage Runge–Kutta time integrators
* [Plots.jl](https://github.com/JuliaPlots/Plots.jl): visualization of the results

They can be added and installed by executing the following commands in the Julia REPL:
```julia
using Pkg
Pkg.add(["TermiteMoundInducedAirflowTrixi", "Trixi", "OrdinaryDiffEqLowStorageRK", "Plots"])
```

## Referencing
You can directly refer to TermiteMoundInducedAirflowTrixi.jl as

```bibtex
@software{marx26,
  title={{TermiteMoundInducedAirflowTrixi.jl}: {T}ermite mound airflow simulations with {T}rixi.jl},
  author={Marx, Oliver P and Gasser, Ingenuin and Schmidgall, Annika},
  year={2026},
  howpublished={\url{https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl}},
  url = {https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl},
  license = {MIT},
  doi = {10.5281/zenodo.22708420}
}
```

When using this package, please make sure to cite **[Trixi.jl](https://github.com/trixi-framework/Trixi.jl)** as well.

## Authors
TermiteMoundInducedAirflowTrixi.jl is maintained by [Oliver P. Marx](https://github.com/oliver-mx).

## License
TermiteMoundInducedAirflowTrixi.jl is licensed under the MIT license (see [LICENSE](LICENSE)).

# TermiteMoundInducedAirflowTrixi

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl/dev/)
[![Build Status](https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl/actions/workflows/CI.yml/badge.svg?)](https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/oliver-mx/TermiteMoundInducedAirflowTrixi.jl/graph/badge.svg?token=w06WNAQpNI)](https://codecov.io/gh/oliver-mx/TermiteMoundInducedAirflowTrixi.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22708420.svg)](https://doi.org/10.5281/zenodo.22708420)

Documentation for [TermiteMoundInducedAirflowTrixi.jl](https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl).

**TermiteMoundInducedAirflowTrixi.jl** is a numerical simulation package written in Julia. It solves a thermo-fluid-dynamic model of the airflow inside a termite mound with the discontinuous Galerkin method. The package builds on the numerical simulation framework for conservation laws [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).

## The termite mound model (TME)

The **1D termite mound model** (termite mound equations, TME) is a thermo-fluid-dynamic model of
the airflow through the flow channel of a termite mound. The channel is a closed loop: it runs
from the soil up through the flute to the top of the mound, down through the chimney and back
through the soil. It is parametrized by the scaled arc length ``x \in [0, 1]`` with periodic
boundary conditions. The model is a low Mach number limit of the Euler equations of gas
dynamics (see [Mathematics](@ref)). All quantities are scaled.

```math
\begin{aligned}
\frac{\partial \rho}{\partial t} + \frac{\partial (\rho u)}{\partial x}
  &= -\frac{A_x}{A} \rho u, \\
\frac{\partial T_i}{\partial t} &= k_i (T - T_i), \\
\frac{\partial p_0}{\partial t} + \gamma p_0 \frac{\partial u}{\partial x}
  &= -\gamma \frac{A_x}{A} u p_0 - \frac{k_w}{\sqrt{A}} (T - T_u), \\
\frac{\mathrm{d} v}{\mathrm{d} t}
  &= \frac{1}{\int_0^1 \frac{\rho}{A} \, dy} \int_0^1 \Big[ -\rho u u_x
     - \frac{\beta}{A} \rho u \big(\eta + (1-\eta)|u|\big) \\
  &\hphantom{= \frac{1}{\int_0^1 \frac{\rho}{A} \, dy} \int_0^1 \Big[}
     - \frac{h_x}{Fr^2} (\rho - \rho_{h_0}) \Big] \, dy, \\
p_0 &= \rho T, \\
u(t, x) &= \frac{v(t)}{A(x)} + \frac{1}{A(x)} \int_0^x \frac{A(y)}{\gamma p_0}
  \left[ -\frac{\mathrm{d} p_0}{\mathrm{d} t} - \frac{k_w}{\sqrt{A}} (T - T_u) \right] dy.
\end{aligned}
```

**Unknowns:**
- air density ``\rho(t, x)`` and air velocity ``u(t, x)``,
- leading-order pressure ``p_0(t)`` (constant in space) and air temperature ``T = p_0/\rho``,
- mound material temperature ``T_i(t, x)``,
- auxiliary velocity ``v(t) = A(0) \, u(t, 0)``, the volume flux at ``x = 0``.

**Model data:**
- cross section ``A(x)`` of the flow channel,
- height profile ``h(x)`` of the flow channel (its slope ``h_x`` drives the buoyancy),
- boundary (wall) temperature ``T_u(t, x, T_i)``: the soil temperature in the soil, the mean
  of the ambient air temperature and ``T_i`` in the flute, and ``T_i`` in the chimney.

The last equation is the **velocity constraint**. Since ``p_0`` does not depend on ``x``, the
``p_0`` equation does not evolve a field. Instead, it determines the velocity ``u``. The model is
therefore a system of a transport PDE for ``\rho``, ODEs for ``v`` and ``p_0``, an ODE for
``T_i`` at every point ``x`` (``T_i`` is not transported) and an algebraic equation for ``u``.

## What the package provides

- **The equations** `TermiteMoundEquations1D` for [Trixi.jl](https://github.com/trixi-framework/Trixi.jl),
  a high-performance framework for conservation laws with the discontinuous Galerkin (DG)
  method, together with the model parameters (`termite_parameters`) and the model data
  (cross section, height profile and boundary temperature) of termite mounds.
- **Two numerical methods** for the PDE-ODE system with algebraic constraint:
  - a **coupled approach** (`CoupledSourceTerms` + `AlgebraicVelocityCallback`): all
    equations are integrated together. **This is the recommended approach.**
  - a **Godunov-type splitting** (`GodunovSplittingCallback`): a first-order operator
    splitting.
- **A manufactured solution** to verify both methods (see [Manufactured solution](@ref)).
- **Example elixirs** in `examples/` for a termite mound, a reference mound and the
  manufactured solution, each for both methods.

To get started, follow the [Tutorial](@ref).

## Installation

If you have not yet installed Julia, please [follow the instructions for your
operating system](https://julialang.org/downloads/platform/). 
TermiteMoundInducedAirflowTrixi.jl works with Julia v1.10 and newer.
We recommend using the latest stable release of Julia.

To run the examples, **TermiteMoundInducedAirflowTrixi.jl**, **Trixi.jl**
(DG solver framework), **OrdinaryDiffEqLowStorageRK.jl** (time integration) and **Plots.jl**
(visualization) can be added and installed by executing the following commands in the Julia REPL:
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

TermiteMoundInducedAirflowTrixi.jl is licensed under the MIT license.

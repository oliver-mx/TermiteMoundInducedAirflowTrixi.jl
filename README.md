# TermiteMoundInducedAirflowTrixi

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://oliver-mx.github.io/TermiteMoundInducedAirflowTrixi.jl/)
[![Build Status](https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl/actions/workflows/CI.yml/badge.svg?)](https://github.com/oliver-mx/TermiteMoundInducedAirflowTrixi.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/oliver-mx/TermiteMoundInducedAirflowTrixi.jl/graph/badge.svg?token=w06WNAQpNI)](https://codecov.io/gh/oliver-mx/TermiteMoundInducedAirflowTrixi.jl)

**TermiteMoundInducedAirflowTrixi.jl** is a numerical simulation package focused on solving a thermo fluid dynamic model describing airflow inside of a termite mound with the discontinuous Galerkin method written in Julia. The package builds on the numerical simulation framework for conservation laws [Trixi.jl](https://github.com/trixi-framework/Trixi.jl). It provides specialized models for termite mound and passive house applications.

## Description

This package provides:

- **1D Termite Mound Model**: This thermo fluid dynamic model describes the airflow dynamics inside of termite mounds. The one-dimensional mathematicalmodel originates from the Euler equations of gas dynamics in a low Mach number regime:

```math
	\begin{equation}
		\frac{\partial \rho}{\partial t} + \frac{\partial \left( \rho u\right)}{\partial x} = - \frac{\textrm{A}_x}{\textrm{A}}\rho u 
	\end{equation}
```
```math
    \begin{equation}
		\frac{\partial v}{\partial t} = \frac{1}{\int_0^1 \frac{\rho}{A} \, dy} \left[ 
		\int_0^1 -\rho u u_x  - \rho u \left( \beta \eta - \beta\left(1-\eta\right) \vert u \vert \right) - \frac{\textrm{h}_x}{Fr^2}(\rho - \rho_{h_0}) \, dy \right] 
	\end{equation}
```
```math
    \begin{equation}
		\frac{\partial p_0}{\partial t} + \gamma p_0 \frac{\partial u}{\partial x}= - \gamma \frac{\textrm{A}_x}{\textrm{A}} u p_0 - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right) 
	\end{equation}
```
```math
    \begin{equation}
		\frac{\partial T_i}{\partial t} = k_i \left( T - T_i \right)
	\end{equation}
```
```math
    \begin{equation}
		p_0 = \rho T 
	\end{equation}
```
```math
    \begin{equation}
		u(t,x) = \frac{v(t)}{A(x)} + \frac{1}{A(x)} \int_0^x A(y)\left[ -\frac{\partial  p_0}{\partial t} - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right)\right] \, dy 
	\end{equation}
```

- **1D Passive House Model**: This model is an extension of the `1D Termite Mound Model`. It also considers radiant flux, therefore modelling the greenhouse effect:
  
```math
	\begin{equation}
		\frac{\partial \rho}{\partial t} + \frac{\partial \left( \rho u\right)}{\partial x} = - \frac{\textrm{A}_x}{\textrm{A}}\rho u 
	\end{equation}
```
```math
    \begin{equation}
		\frac{\partial v}{\partial t} = \frac{1}{\int_0^1 \frac{\rho}{A} \, dy} \left[ 
		\int_0^1 -\rho u u_x  - \rho u \left( \beta \eta - \beta\left(1-\eta\right) \vert u \vert \right) - \frac{\textrm{h}_x}{Fr^2}(\rho - \rho_{h_0}) \, dy \right] 
	\end{equation}
```
```math
    \begin{equation}
		\frac{\partial p_0}{\partial t} + \gamma p_0 \frac{\partial u}{\partial x}= - \gamma \frac{\textrm{A}_x}{\textrm{A}} u p_0 - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right) 
	\end{equation}
```
```math
    \begin{equation}
		\frac{\partial T_i}{\partial t} = k_i \left( T - T_i \right)  + k_s \textrm{q}_{\textrm{s}} 
	\end{equation}
```
```math
    \begin{equation}
		p_0 = \rho T 
	\end{equation}
```
```math
    \begin{equation}
		u(t,x) = \frac{v(t)}{A(x)} + \frac{1}{A(x)} \int_0^x A(y)\left[ -\frac{\partial  p_0}{\partial t} - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right)\right] \, dy 
	\end{equation}
```

More details about the models can be found in the publication:  **[[Article](...)]**

Both models are designed to be used with **[Trixi.jl](https://github.com/trixi-framework/Trixi.jl)**, a flexible and high-performance framework for solving systems of conservation laws using the Discontinuous Galerkin (DG) method.

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

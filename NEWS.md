# Changelog

TermiteMoundInducedAirflowTrixi.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file for human readability.

## Changes when updating to v2.0.0 from v1.1.x

#### Added

- "Coupled" method: `CoupledSourceTerms` and `AlgebraicVelocityCallback`.
- Model data: `CrossSectionConstant`, `CrossSectionSine`, `HeightProfileTermiteMound`,
  `BoundaryTemperatureTermiteMound`.
- Manufactured solution: `ManufacturedInitialCondition`, `manufactured_source_terms`,
  `manufactured_source_terms_ode`.
- New examples and documentation pages.

#### Changed

- `UpdateVelocityCallback` was renamed to `GodunovSplittingCallback`.
- `TermiteMoundEquations1D` requires `cross_section`, `height_profile` and
  `boundary_temperature`.
- FastGaussQuadrature.jl is no longer a dependency.

#### Fixed

- Sign in the friction term: `η + (1 - η)|u|` instead of `η - (1 - η)|u|`.

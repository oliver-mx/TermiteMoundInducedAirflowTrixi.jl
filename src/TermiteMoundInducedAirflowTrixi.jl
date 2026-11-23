module TermiteMoundInducedAirflowTrixi

using Trixi
using OrdinaryDiffEqLowStorageRK
using Interpolations
# Import additional symbols that are not exported by Trixi.jl
using MuladdMacro: @muladd
using Trixi: AbstractEquations
import Interpolations: Line
# (all other extended Trixi functions are written fully qualified, e.g. `Trixi.flux`)
import Trixi: flux_ranocha, varnames, ln_mean

include("equations/equations.jl")
include("solvers/solvers.jl")
include("callbacks_step/callbacks_step.jl")

# Export types/functions that define the public API of TermiteMoundInducedAirflowTrixi.jl
export TermiteMoundEquations1D

export TermiteMoundInitialCondition, source_terms, termite_parameters

# model data: cross section, height profile and boundary temperature
export CrossSectionConstant, CrossSectionSine, CrossSection
export HeightProfileTermiteMound, HeightProfile
export BoundaryTemperatureTermiteMound, BoundaryTemperature

# "splitting" method:
export GodunovSplittingCallback

# "coupled" method: coupled method of lines
export CoupledSourceTerms, AlgebraicVelocityCallback

export ManufacturedInitialCondition,
    manufactured_source_terms, manufactured_source_terms_ode

end

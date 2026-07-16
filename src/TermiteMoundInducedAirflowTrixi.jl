module TermiteMoundInducedAirflowTrixi

using Trixi
using OrdinaryDiffEqLowStorageRK
using Interpolations
using QuadGK
using FastGaussQuadrature
# Import additional symbols that are not exported by Trixi.jl
using MuladdMacro: @muladd
using Trixi: AbstractEquations
import Interpolations: Line
import Trixi: flux_ranocha, ln_mean, inv_ln_mean, flux, varnames, cons2cons, cons2prim, prim2cons, cons2entropy, max_abs_speeds

include("equations/equations.jl")
include("callback_step/callback_step.jl")

# Export types/functions that define the public API of TermiteMoundInducedAirflowTrixi.jl
export TermiteMoundEquations1D, PassiveHouseEquations1D

export UpdateVelocityCallback

end
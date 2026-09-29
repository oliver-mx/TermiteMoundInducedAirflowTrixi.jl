module TermiteMoundInducedAirflowTrixi

using Trixi
using OrdinaryDiffEqLowStorageRK
using Interpolations
using FastGaussQuadrature
# Import additional symbols that are not exported by Trixi.jl
using MuladdMacro: @muladd
using Trixi: AbstractEquations
import Interpolations: Line
import Trixi:
    flux_ranocha,
    flux_hllc,
    ln_mean,
    inv_ln_mean,
    flux,
    varnames,
    cons2prim,
    prim2cons,
    cons2entropy,
    cons2cons,
    max_abs_speeds,
    max_abs_speed_naive

include("equations/equations.jl")
include("callback_step/callback_step.jl")

# Export types/functions that define the public API of TermiteMoundInducedAirflowTrixi.jl
export TermiteMoundEquations1D

export UpdateVelocityCallback,
    TermiteMoundInitialCondition, source_terms, termite_parameters

end

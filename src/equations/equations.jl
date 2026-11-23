# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

####################################################################################################
# Include files with actual implementations for different systems of equations. 

# model data passed to the equations: cross section, height profile, boundary temperature
include("model_data.jl")
# Termite mound equations (TME) in 1D: equations, parameters, initial condition, source terms,
# numerical fluxes
include("termite_mound_1d.jl")
# manufactured solution of the forced TME: exact fields, source terms of all equations
include("manufactured_solution.jl")
end # @muladd

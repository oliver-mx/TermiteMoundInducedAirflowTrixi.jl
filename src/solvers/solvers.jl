# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

####################################################################################################
# Include files with additional solver functionality.

# DG-consistent nodal operators, shared by the "splitting" and the "coupled" methods
include("dg_operators.jl")
# right-hand side of the "coupled" method
include("coupled_rhs.jl")
end # @muladd

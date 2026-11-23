# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

####################################################################################################
# Include files with callbacks that are called after every time step.

# "splitting" method: first-order Godunov operator splitting
include("godunov_splitting.jl")
# "coupled" method: algebraic velocity for step size control, analysis and output
include("algebraic_velocity.jl")
end # @muladd

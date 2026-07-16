using TermiteMoundInducedAirflowTrixi
using Trixi, Trixi2Vtk, OrdinaryDiffEqLowStorageRK, Interpolations, QuadGK, FastGaussQuadrature, Plots
using Trixi: AbstractEquations, @muladd
import Interpolations: Line
import Trixi: flux_ranocha, ln_mean, inv_ln_mean, flux, varnames, cons2cons, cons2prim, prim2cons, cons2entropy, max_abs_speeds
import TermiteMoundInducedAirflowTrixi: TermiteMoundInitialCondition, source_terms

###############################################################################
# Semidiscretization

equations = TermiteMoundEquations1D()

initial_condition = TermiteMoundInitialCondition
volume_flux = flux_ranocha
surface_flux = flux_ranocha

dg = DGSEM(polydeg = 3, surface_flux = flux_ranocha, volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

mesh = TreeMesh((0.0,), (1.0,), initial_refinement_level = 5, n_cells_max = 1000, periodicity = true) 

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg, source_terms = source_terms);

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
stepsize_callback = StepsizeCallback(cfl = 0.5)
update_velocity_callback = UpdateVelocityCallback()

callbacks = CallbackSet(summary_callback, stepsize_callback, update_velocity_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false); dt = 1.0,
            ode_default_options()..., callback = callbacks);

#pd = PlotData1D((x, equations) -> initial_condition(x, last(tspan), equations), semi); plot(pd)
#pd = PlotData1D(sol); plot(pd);

using TermiteMoundInducedAirflowTrixi
using Trixi, OrdinaryDiffEqLowStorageRK

###############################################################################
# Semidiscretization
#
# "coupled" method (coupled method of lines, `CoupledSourceTerms`): all equations are integrated
# together; the velocity is computed from the constraint in every RHS evaluation.

(γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr², ρₕ₀, T_ref, t_ref, T0, v0, Ti_LI) =
    termite_parameters(0.6, 2.0)
# model data: constant cross section A = 1, termite mound height profile and boundary
# temperature
cross_section = CrossSectionConstant()
height_profile = HeightProfileTermiteMound(; xa, xb, L, h)
boundary_temperature = BoundaryTemperatureTermiteMound(; xa, xb, xc, tᵣ, T_ref)
equations = TermiteMoundEquations1D(;
    γ,
    k_i,
    k_w,
    tᵣ,
    uᵣ,
    xa,
    xb,
    xc,
    r,
    h,
    L,
    β,
    η,
    Fr²,
    ρₕ₀,
    T_ref,
    t_ref,
    T0,
    v0,
    Ti_LI,
    cross_section,
    height_profile,
    boundary_temperature,
)

initial_condition = TermiteMoundInitialCondition
volume_flux = flux_ranocha
surface_flux = flux_hllc

dg = DGSEM(
    polydeg = 4,
    surface_flux = surface_flux,
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux),
)

mesh = TreeMesh(
    (0.0,),
    (1.0,),
    initial_refinement_level = 5,
    n_cells_max = 1000,
    periodicity = true,
)

semi = SemidiscretizationHyperbolic(
    mesh,
    equations,
    initial_condition,
    dg,
    source_terms = CoupledSourceTerms(source_terms),
    boundary_conditions = boundary_condition_periodic,
);

###############################################################################
# ODE solvers, callbacks etc.

# one full day; the analysis "errors" compare the final state with the flat initial condition
tspan = (0.0, 1.0) .* 86400 ./ equations.tᵣ
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
algebraic_velocity_callback = AlgebraicVelocityCallback()
stepsize_callback = StepsizeCallback(cfl = 3.0)
analysis_callback = AnalysisCallback(semi, interval = 100_000, uEltype = real(dg))

callbacks = CallbackSet(
    summary_callback,
    algebraic_velocity_callback,
    stepsize_callback,
    analysis_callback,
)

###############################################################################
# run the simulation

sol = solve(
    ode,
    CarpenterKennedy2N54(williamson_condition = false);
    dt = 1.0,
    ode_default_options()...,
    callback = callbacks,
);

#────────────────────────────────────────────────────────────────────────────────────────────────────
# Simulation running 'TermiteMoundEquations1D' with DGSEM(polydeg=4) and cfl=0.05
#────────────────────────────────────────────────────────────────────────────────────────────────────
#timesteps:            4981260                run time:       5.38328991e+02 s
# Δt:             1.00243664e-04                └── GC time:    5.10689598e+01 s (9.487%)
# sim. time:      8.53943670e+02 (100.000%)     time/DOF/rhs!:  2.15988670e-08 s
#                                               PID:            1.34832881e-07 s
#DOFs per field:           160                alloc'd memory:         61.473 MiB
#elements:                  32                device memory:           0.000 MiB
#
# Variable:       rho              v1               p                Ti
# L2 error:       2.03338087e-03   2.64640909e-01   1.95860867e-03   2.73873479e-04
# Linf error:     3.97670479e-03   2.67574474e-01   1.95860867e-03   4.67148540e-04
# ∑∂S/∂U ⋅ Uₜ :  -9.54230560e-07
#────────────────────────────────────────────────────────────────────────────────────────────────────

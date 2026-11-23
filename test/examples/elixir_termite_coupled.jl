using TermiteMoundInducedAirflowTrixi
using Trixi, OrdinaryDiffEqLowStorageRK

###############################################################################
# Semidiscretization
#
# "coupled" method (coupled method of lines, `CoupledSourceTerms`): all equations are integrated
# together; the velocity is computed from the constraint in every RHS evaluation.

(γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr², ρₕ₀, T_ref, t_ref, T0, v0, Ti_LI) =
    termite_parameters(1.5, 2.5)
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
#timesteps:            3384561                run time:       3.86713752e+02 s
# Δt:             9.96086634e-05                └── GC time:    3.82645674e+01 s (9.895%)
# sim. time:      5.96258330e+02 (100.000%)     time/DOF/rhs!:  2.26147231e-08 s
#                                               PID:            1.44054609e-07 s
#DOFs per field:           160                alloc'd memory:         67.491 MiB
#elements:                  32                device memory:           0.000 MiB
#
# Variable:       rho              v1               p                Ti
# L2 error:       2.76070708e-03   3.10680275e-01   2.28637208e-03   8.99294108e-05
# Linf error:     4.73205376e-03   3.13950078e-01   2.28637208e-03   1.88372920e-04
# ∑∂S/∂U ⋅ Uₜ :  -9.38332591e-08
#────────────────────────────────────────────────────────────────────────────────────────────────────

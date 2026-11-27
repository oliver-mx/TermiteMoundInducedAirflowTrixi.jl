using TermiteMoundInducedAirflowTrixi
using Trixi, OrdinaryDiffEqLowStorageRK

###############################################################################
# Semidiscretization
#
# "splitting" method: first-order Godunov (Lie-Trotter) operator splitting (`GodunovSplittingCallback`).

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
    source_terms = source_terms,
    boundary_conditions = boundary_condition_periodic,
);

###############################################################################
# ODE solvers, callbacks etc.

# one full day; the analysis "errors" compare the final state with the flat initial condition
tspan = (0.0, 1.0) .* 86400 ./ equations.tᵣ
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
stepsize_callback = StepsizeCallback(cfl = 3.0)
godunov_splitting_callback =
    GodunovSplittingCallback(CarpenterKennedy2N54(williamson_condition = false))
analysis_callback = AnalysisCallback(semi, interval = 100_000, uEltype = real(dg))

callbacks = CallbackSet(
    summary_callback,
    stepsize_callback,
    godunov_splitting_callback,
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
#timesteps:            3423744                run time:       3.74372513e+02 s
# Δt:             1.02120765e-04                └── GC time:    3.63211302e+01 s (9.702%)
# sim. time:      5.96258330e+02 (100.000%)     time/DOF/rhs!:  2.15026301e-08 s
#                                               PID:            1.35707484e-07 s
#DOFs per field:           160                alloc'd memory:         76.465 MiB
#elements:                  32                device memory:           0.000 MiB
#
# Variable:       rho              v1               p                Ti
# L2 error:       2.85065479e-03   2.62533517e-01   2.46944669e-03   9.18868740e-05
# Linf error:     4.81726039e-03   2.66092963e-01   2.46944669e-03   1.76871225e-04
# ∑∂S/∂U ⋅ Uₜ :  -4.78897747e-08
#────────────────────────────────────────────────────────────────────────────────────────────────────

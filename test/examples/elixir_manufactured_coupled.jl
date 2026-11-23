using TermiteMoundInducedAirflowTrixi
using Trixi, OrdinaryDiffEqLowStorageRK

###############################################################################
# Semidiscretization
#
# Manufactured solution with variable cross section A(x) = 1 + 0.3 sin(2πx + π/2) (see
# `ManufacturedInitialCondition`, `A`): the analysis errors are the
# deviation of the numerical solution from the exact solution. With the "coupled" method (coupled
# method of lines, `CoupledSourceTerms`), the manufactured sources of rho and Ti (`manufactured_source_terms`)
# and of p0 and v (`manufactured_source_terms_ode`) enter one coupled RHS.

(γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr², ρₕ₀, T_ref, t_ref, T0, v0, Ti_LI) =
    termite_parameters(1.5, 2.5)
# model data: variable cross section A(x) = 1 + 0.3 sin(2πx + π/2), termite mound height
# profile and boundary temperature
cross_section = CrossSectionSine(0.3)
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

initial_condition = ManufacturedInitialCondition
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
    source_terms = CoupledSourceTerms(
        manufactured_source_terms;
        ode_source_terms = manufactured_source_terms_ode,
    ),
    boundary_conditions = boundary_condition_periodic,
);

###############################################################################
# ODE solvers, callbacks etc.

# one time period of the manufactured solution
tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
algebraic_velocity_callback = AlgebraicVelocityCallback()
stepsize_callback = StepsizeCallback(cfl = 1.0)
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
